/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "../offload/offload_library.h"
#include "../offload/offload_runtime.h"
#include "dbm_mempool.h"
#include "dbm_mpi.h"

/*******************************************************************************
 * \brief Private routine for actually allocating system memory.
 * \author Ole Schuett
 ******************************************************************************/
static void *actual_malloc(const size_t size, const bool on_device) {
  void *memory = NULL;

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  if (on_device) {
    offload_activate_chosen_device();
    offloadMalloc(&memory, size);
  } else {
    offloadMallocHost(&memory, size);
  }
#else
  (void)on_device; // mark used
  memory = dbm_mpi_alloc_mem(size);
#endif

  assert(memory != NULL);
  return memory;
}

/*******************************************************************************
 * \brief Private routine for actually freeing system memory.
 * \author Ole Schuett
 ******************************************************************************/
static void actual_free(void *memory, const bool on_device) {
  if (memory == NULL) {
    return;
  }

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  if (on_device) {
    offload_activate_chosen_device();
    offloadFree(memory);
  } else {
    offloadFreeHost(memory);
  }
#else
  (void)on_device; // mark used
  dbm_mpi_free_mem(memory);
#endif
}

/*******************************************************************************
 * \brief Private struct for storing a chunk of memory.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct dbm_memchunk {
  void *mem; // first: allows to cast memchunk into mem-ptr...
  struct dbm_memchunk *next;
  size_t size;
  bool on_device;
} dbm_memchunk_t;

/*******************************************************************************
 * \brief Private linked list of memory chunks that are available.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_memchunk_t *mempool_available_head = NULL;

/*******************************************************************************
 * \brief Private linked list of memory chunks that are in use.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_memchunk_t *mempool_allocated_head = NULL;

/*******************************************************************************
 * \brief Private statistics (survives dbm_mempool_clear).
 * \author Hans Pabst
 ******************************************************************************/
static dbm_memstats_t mempool_stats = {0};

/*******************************************************************************
 * \brief Private routine for allocating host or device memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
static void *internal_mempool_malloc(const size_t size, const bool on_device) {
  if (size == 0) {
    return NULL;
  }

  dbm_memchunk_t *chunk = NULL;
#pragma omp critical(dbm_mempool_modify)
  {
    // Find a suitable chunk in mempool_available.
    dbm_memchunk_t **indirect = &mempool_available_head;
    dbm_memchunk_t **hit = NULL, **fallback = NULL;
    for (; NULL != *indirect; indirect = &(*indirect)->next) {
      if ((*indirect)->on_device == on_device) {
        if ((*indirect)->size < size) {
          if (NULL == fallback || (*fallback)->size < (*indirect)->size) {
            fallback = indirect;
          }
        } else if (NULL == hit || (*indirect)->size < (*hit)->size) {
          hit = indirect;
          if (size == (*hit)->size) {
            break;
          }
        }
      }
    }
    if (NULL == hit) {
      hit = fallback;
    }

    // If a chunck was found, remove it from mempool_available.
    if (NULL != hit) {
      chunk = *hit;
      *hit = chunk->next;
      assert(chunk->on_device == on_device);
    } else { // Allocate a new chunk.
      assert(chunk == NULL);
      chunk = malloc(sizeof(dbm_memchunk_t));
      assert(chunk != NULL);
      chunk->on_device = on_device;
      chunk->size = 0;
      chunk->mem = NULL;
    }

    // Insert chunk into mempool_allocated.
    chunk->next = mempool_allocated_head;
    mempool_allocated_head = chunk;

    if (chunk->size < size) {
      // Update statistics before resizing chunk
      if (on_device) {
        mempool_stats.device_size += size - chunk->size;
        ++mempool_stats.device_mallocs;
      } else {
        mempool_stats.host_size += size - chunk->size;
        ++mempool_stats.host_mallocs;
      }
      // Resize chunk if needed
      actual_free(chunk->mem, chunk->on_device);
      chunk->mem = actual_malloc(size, chunk->on_device);
      chunk->size = size; // update
    }
  }

  return chunk->mem;
}

/*******************************************************************************
 * \brief Internal routine for allocating host memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
void *dbm_mempool_host_malloc(const size_t size) {
  return internal_mempool_malloc(size, false);
}

/*******************************************************************************
 * \brief Internal routine for allocating device memory from the pool
 * \author Ole Schuett
 ******************************************************************************/
void *dbm_mempool_device_malloc(const size_t size) {
  return internal_mempool_malloc(size, true);
}

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_free(void *mem) {
  if (mem == NULL) {
    return;
  }

#pragma omp critical(dbm_mempool_modify)
  {
    // Find chunk in mempool_allocated.
    dbm_memchunk_t **indirect = &mempool_allocated_head;
    while (*indirect != NULL && (*indirect)->mem != mem) {
      indirect = &(*indirect)->next;
    }
    dbm_memchunk_t *chunk = *indirect;
    assert(chunk != NULL && chunk->mem == mem);

    // Remove chunk from mempool_allocated.
    *indirect = chunk->next;

    // Add chunk to mempool_available.
    chunk->next = mempool_available_head;
    mempool_available_head = chunk;
  }
}

/*******************************************************************************
 * \brief Internal routine for freeing all memory in the pool (not thread-safe).
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_clear(void) {
  assert(omp_get_num_threads() == 1);
  assert(mempool_allocated_head == NULL); // check for memory leak
#if 0
  while (mempool_allocated_head != NULL) {
    dbm_memchunk_t *chunk = mempool_allocated_head;
    mempool_allocated_head = chunk->next;
    printf("Found alloacted memory chunk of size: %lu\n", chunk->size);
    actual_free(chunk->mem, chunk->on_device);
    free(chunk);
  }
#endif
  // Free chunks in mempool_available.
  while (mempool_available_head != NULL) {
    dbm_memchunk_t *chunk = mempool_available_head;
    mempool_available_head = chunk->next;
    actual_free(chunk->mem, chunk->on_device);
    free(chunk);
  }
}

/*******************************************************************************
 * \brief Internal routine to query statistics (not thread-safe).
 * \author Hans Pabst
 ******************************************************************************/
void dbm_mempool_statistics(dbm_memstats_t *memstats) {
  assert(NULL != memstats);
  *memstats = mempool_stats;
}

// EOF
