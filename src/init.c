#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void R_init_gremlin(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}

