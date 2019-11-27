#include <stdio.h>

#ifndef BUFFERS_C_
#define BUFFERS_C_

// ================================================
// Buffer for keeping track of records
// ================================================
typedef struct {
    uint32_t *records;
    size_t head;
} record_buf_t;


void record_buf_reset(record_buf_t *buffer) {
    buffer->head = 0;
}

// ================================================
// END Buffer for keeping track of records
// ================================================

// ====================================================================
// ring BUFFER Structure for use with the g2 algorithm based on it
// ====================================================================
typedef struct {
    uint64_t *buffer;
    int head;
    int count;
    int size; //of the buffer
} ring_buf_t;


void ring_buf_reset(ring_buf_t * cbuf);
static inline void ring_buf_put(ring_buf_t * cbuf, uint64_t data);
static inline uint64_t ring_buf_oldest(ring_buf_t * cbuf);

ring_buf_t ring_buf_allocate(int size) {
    ring_buf_t cbuf;
    cbuf.size = size;
    cbuf.buffer = malloc(cbuf.size * sizeof(uint64_t)); // set memory to zero so we have a proper
    ring_buf_reset(&cbuf);

    return cbuf;
}

void ring_buf_reset(ring_buf_t * cbuf)
{
    if(cbuf) {
        cbuf->head = 0;
        cbuf->count = 0;
        for (int i = 0; i < cbuf->size; ++i)
        {
            cbuf->buffer[i] = 0;
        }
    }
}

static inline void ring_buf_put(ring_buf_t * cbuf, uint64_t data)
{
    cbuf->buffer[cbuf->head] = data;
    cbuf->head = (cbuf->head + 1) % cbuf->size;
    if(cbuf->count < cbuf->size) {
        cbuf->count = cbuf->count + 1;
    }
}

static inline uint64_t ring_buf_oldest(ring_buf_t * cbuf) {
    // CAUTION: Even if the buffer is empty oldest will return whatever
    // is in buffer[0]. We do so because it is conveninet for our specific
    // application but it can be catastrophic.
    // TAKE HOME MESSAGE: Don't use this implementation as is for anything
    // other than computing g2.
    if(cbuf->count < cbuf->size) {
        return cbuf->buffer[0];
    }
    return cbuf->buffer[cbuf->head];
}

static inline void ring_buf_grow(ring_buf_t *cbuf) {
//    if (cbuf->size < 2048*2048*2048)
//    {   
        cbuf->size = 2*cbuf->size;  // double the space as previous array
        cbuf->buffer = realloc(cbuf->buffer, cbuf->size * sizeof(uint64_t));
        // set memory to zero so we have a proper
        for (int i = (cbuf->size)/2; i < cbuf->size; ++i)
        {
            cbuf->buffer[i]=0;
        }
//    }
}

// ============================
// END OF ring BUFFER
// ============================

// ====================================================================
// CHAINED LIST BUFFER Structure for use with the dummy g2 algorithm (calculate_g2)
// ====================================================================

// this 'node' structure is a chained list to be used with the buffers. It is easy to add an item
// at the end and read+remove an item at the beginning (see the three functions below).
typedef struct node {
    uint64_t val;
    struct node * next;
} node_t;

void head_init(node_t** head, int* length) {
    // initialize a chained list with value 0 and length 0.
    *head = malloc(sizeof(node_t));
    (*head)->val = 0;
    (*head)->next = NULL;
    *length = 0;
}

void push(node_t * head, uint64_t val, int* length) {
    // add an item at the end of the list
    node_t * current = head;
    while (current->next != NULL) {
        current = current->next;
    }
    
    // now we can add a new variable
    current->next = malloc(sizeof(node_t));
    current->next->val = val;
    current->next->next = NULL;
    *length += 1;
}

uint64_t pop(node_t * head, int* length) {
    // remove the first info item (second in the list) from the list, returning its value
    // remember that, for simplicity, we keep the very first item of the list
    // alive so we don't need to recreate one each time the list becomes empty.
    uint64_t retval = 0;
    node_t * next_node = NULL;
    
    if (head->next == NULL) {
        return 0;
    }
    
    next_node = head->next->next;
    retval = head->next->val;
    free(head->next);
    head->next = next_node;
    *length = *length - 1;
    
    return retval;
}

// ============================
// END OF CHAINED LIST BUFFER
// ============================

#endif /* BUFFERS_C_ */