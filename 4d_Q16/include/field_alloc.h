// -----------------------------------------------------------------
// Macros to allocate space for fields
#ifndef _FIELD_ALLOC_H
#define _FIELD_ALLOC_H

#define FIELD_ALLOC(name, typ) {                        \
  name = malloc(sites_on_node * sizeof(typ));           \
  if (name == NULL) {                                   \
    printf("node%d: FIELD_ALLOC failed\n", this_node);  \
    terminate(1);                                       \
  }                                                     \
}

#define FIELD_ALLOC_VEC(name, typ, size) {                    \
  int ifield;                                                 \
  for (ifield = 0; ifield < size; ifield++) {                 \
    name[ifield] = malloc(sites_on_node * sizeof(typ));       \
    if (name[ifield] == NULL) {                               \
      printf("node%d: FIELD_ALLOC_VEC failed\n", this_node);  \
      terminate(1);                                           \
    }                                                         \
  }                                                           \
}

#define FIELD_ALLOC_MAT(name, typ, sizea, sizeb) {                      \
  int ifield, jfield;                                                   \
  for (ifield = 0; ifield < sizea; ifield++) {                          \
    for (jfield = 0; jfield < sizeb; jfield++) {                        \
      name[ifield][jfield] = malloc(sites_on_node * sizeof(typ));       \
      if (name[ifield][jfield] == NULL) {                               \
        printf("node%d: FIELD_ALLOC_MAT_OFFDIAG failed\n", this_node);  \
        terminate(1);                                                   \
      }                                                                 \
    }                                                                   \
  }                                                                     \
}

#define FIELD_ALLOC_MAT_OFFDIAG(name, typ, size) {                        \
  int ifield, jfield;                                                     \
  for (ifield = 0; ifield < size; ifield++) {                             \
    for (jfield = 0; jfield < size; jfield++) {                           \
      if (ifield != jfield) {                                             \
        name[ifield][jfield] = malloc(sites_on_node * sizeof(typ));       \
        if (name[ifield][jfield] == NULL) {                               \
          printf("node%d: FIELD_ALLOC_MAT_OFFDIAG failed\n", this_node);  \
          terminate(1);                                                   \
        }                                                                 \
      }                                                                   \
    }                                                                     \
  }                                                                       \
}
#endif
// -----------------------------------------------------------------
