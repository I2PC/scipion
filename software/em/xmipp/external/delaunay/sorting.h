#ifndef SORTING_H
#define SORTING_H

#include "dcel.h"


/****************************************************************************
                           PUBLIC FUNCTION HEADERS
*****************************************************************************/
void	sort(struct DCEL_T *dcel);
void	clutter(struct DCEL_T *dcel);
void    set_Highest_First(struct DCEL_T *dcel, int (*f)(struct Point_T *, struct Point_T *));

#endif
