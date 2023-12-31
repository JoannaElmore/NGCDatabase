/* ==================================================== */
/*      Copyright (c) Atelier de BioInformatique        */
/*      Mar. 92                                         */
/*      File: libstki.h                                 */
/*      Purpose: library of dynamic stacks holding      */
/*               integer values                         */
/*      History:                                        */
/*      00/03/92 : <Gloup> first draft                  */
/*      07/07/93 : <Gloup> complete revision            */
/*      10/03/94 : <Gloup> added xxxVector funcs        */
/*      14/05/99 : <Gloup> last revision                */
/* ==================================================== */

#ifndef _H_Gtypes
#include "Gtypes.h"
#endif

#define _H_libstki

/* ==================================================== */
/* Constantes de dimensionnement                        */
/* ==================================================== */

#ifndef kMinStackiSize
#define kMinStackiSize          2       /* taille mini stack    */
#endif


#define kStkiNoErr              0       /* ok                   */
#define kStkiMemErr             1       /* not enough memory    */

#define kStkiReset              Vrai
#define kStkiGet                Faux

/* ==================================================== */
/* Macros standards                                     */
/* ==================================================== */

#ifndef NEW
#define NEW(typ)                (typ*)malloc(sizeof(typ)) 
#define NEWN(typ, dim)          (typ*)malloc((unsigned long)(dim) * sizeof(typ))
#define REALLOC(typ, ptr, dim)  (typ*)realloc((void *) (ptr), (unsigned long)(dim) * sizeof(typ))
#define FREE(ptr)               free((Ptr) ptr)
#endif


/* ==================================================== */
/*  Types & Structures de donnees                       */
/* ==================================================== */

                                                /* -------------------- */
                                                /* structure : pile     */
                                                /* -------------------- */
typedef struct Stacki {
                                                /* ---------------------*/
        Int32                   size;           /* stack size           */
        Int32                   top;            /* current free pos.    */
        Int32                   cursor;         /* current cursor       */
        Int32                   *val;           /* values               */
                                                /* ---------------------*/
} Stacki, *StackiPtr, **StackiHdle;      



/* ==================================================== */
/*  Prototypes (generated by mproto)                    */
/* ==================================================== */

                                        /* libstki.c    */
                                        
Int16           StkiError          (Bool       reset                    );
StackiPtr       NewStacki          (Int32      size                     );
StackiPtr       FreeStacki         (StackiPtr  stki                     );
StackiHdle      NewStackiVector    (Int32 vectSize,   Int32 stackSize   );
StackiHdle      FreeStackiVector   (StackiHdle stkh,  Int32 vectSize    );
Int32           ResizeStacki       (StackiHdle stkh , Int32 size        );
Bool            PushiIn            (StackiHdle stkh , Int32 val         );
Bool            PopiOut            (StackiHdle stkh , Int32 *val        );
Bool            ReadiDown          (StackiPtr  stki , Int32 *val        );
Bool            ReadiUp            (StackiPtr  stki , Int32 *val        );
void            CursiToTop         (StackiPtr  stki                     );
void            CursiToBottom      (StackiPtr  stki                     );
void            CursiSwap          (StackiPtr  stki                     );
Bool            SearchDownStacki   (StackiPtr  stki  , Int32 sval       );
Bool            BinSearchStacki    (StackiPtr  stki  , Int32 sval       );
Bool            SameStacki         (StackiPtr  stki1 , StackiPtr stki2  );
Bool            ReverseStacki      (StackiPtr  stki                     );
