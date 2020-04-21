#ifndef SQL_CALLBACK_HPP
#define SQL_CALLBACK_HPP

#include <sqlite3.h> 

struct dbEntry {
  std::string id = "none";
  float k   = 0.0;
  float g   = 0.0;
  float s   = 0.0;
  float dkg = 0.0;
  float ds  = 0.0;
};

static int callback(void* userData,int argc,char** argv,char** azColName){
  if( argc > 0 ){
    dbEntry* E = (dbEntry*) userData;
    for(int i=0;i<argc;i++){
      if( strcmp(azColName[i],"id") == 0 ){
	E->id = argv[i];
      } else if( strcmp(azColName[i],"k") == 0 ){
	E->k = atof(argv[i]);
      } else if( strcmp(azColName[i],"g") == 0 ){
	E->g = atof(argv[i]);
      } else if( strcmp(azColName[i],"s") == 0 ){
	E->s = atof(argv[i]);
      } else if( strcmp(azColName[i],"dkg") == 0 ){
	E->dkg = atof(argv[i]);
      } else {
	E->ds = atof(argv[i]);
      }
    }
  }
  return 0;
}

#endif /* SQL_CALLBACK_HPP */
