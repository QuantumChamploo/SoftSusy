/* Functions needed for fermion self energies */

#ifndef _SELF_FERMION_H_
#define _SELF_FERMION_H_

/* #include "internal.h" */

int bFS (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL,
         TSIL_COMPLEX *, TSIL_COMPLEX *);

int bpFS (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL,
          TSIL_COMPLEX *, TSIL_COMPLEX *);

int bFV (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL,
         TSIL_COMPLEX *, TSIL_COMPLEX *);

int mSFFSF (TSIL_DATA *, int, 
            TSIL_COMPLEX *, TSIL_COMPLEX *,
            TSIL_COMPLEX *, TSIL_COMPLEX *,
            TSIL_COMPLEX *, TSIL_COMPLEX *);

int mSSFFS (TSIL_DATA *, 
            TSIL_COMPLEX *, TSIL_COMPLEX *, TSIL_COMPLEX *);

int yFSSS (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL,
           TSIL_REAL, TSIL_REAL,
           TSIL_COMPLEX *, TSIL_COMPLEX *);

int vFSSSS (TSIL_DATA *, TSIL_DATA *,
            TSIL_COMPLEX *, TSIL_COMPLEX *);

int vSFFFS (TSIL_DATA *, TSIL_DATA *, int,
            TSIL_COMPLEX *, TSIL_COMPLEX *,
            TSIL_COMPLEX *, TSIL_COMPLEX *,
            TSIL_COMPLEX *, TSIL_COMPLEX *);

int vFSSFF (TSIL_DATA *, TSIL_DATA *,
            TSIL_COMPLEX *, TSIL_COMPLEX *, 
            TSIL_COMPLEX *, TSIL_COMPLEX *);

int f1f4 (TSIL_REAL, TSIL_REAL, TSIL_REAL z, TSIL_REAL, 
          TSIL_COMPLEX *, TSIL_COMPLEX *);

int f2f3f5f6 (TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, 
              TSIL_COMPLEX *, TSIL_COMPLEX *,
              TSIL_COMPLEX *, TSIL_COMPLEX *);

TSIL_REAL F1 (TSIL_REAL, TSIL_REAL);

TSIL_REAL F2 (TSIL_REAL, TSIL_REAL);

TSIL_REAL F3 (TSIL_REAL, TSIL_REAL, TSIL_REAL);

TSIL_REAL F4 (TSIL_REAL, TSIL_REAL, TSIL_REAL);

TSIL_REAL f (TSIL_REAL);

#endif /* self_fermion.h */