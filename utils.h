#ifndef _UTILS_H_
#define _UTILS_H_

#include <sqlite3.h>

/* step when only one row is expected */
void sqlite3_step_onerow(sqlite3_stmt *stmt);

#endif /* _UTILS_H_ */
