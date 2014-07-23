#include <stdlib.h>
#include <stdio.h>

#include "utils.h"
#include "errors.h"

void sqlite3_step_onerow(sqlite3_stmt *stmt)
{
     int res_code;
     res_code = sqlite3_step(stmt);
     if (res_code == SQLITE_DONE) {
          fprintf(stderr, "No rows returned for query: %s\n",
                  sqlite3_sql(stmt));
     }
     if (res_code != SQLITE_ROW) {
          fprintf(stderr, "SQLite3 error: %d (%s)\n",
                  res_code, sqlite3_sql(stmt));
          exit(SQL_ERROR);
     }
}
