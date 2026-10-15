#include <fstream>
#include <sqlite3.h> 
#include <string>
#include <cmath>
#include "json/json.h"

#include "sql_callback.hpp"

int main(int argc, char* argv[]) {

  std::string dbfile = argv[1];
  std::string out_path = argv[2];
  
  // Read the multiple images' parameters from JSON
  Json::Value images;
  std::ifstream fin(out_path+"output/multiple_images.json",std::ifstream::in);
  fin >> images;
  fin.close();


  double kg_sep = 0.05;

  
  double dtmax = 0;
  int i_dtmax = 0;
  for(int q=0;q<images.size();q++){
    if( images[q]["dt"].asDouble() > dtmax ){
      dtmax = images[q]["dt"].asDouble();
      i_dtmax = q;
    }
  }

    
  sqlite3* db;
  char* zErrMsg = 0;
  int rc;
  rc = sqlite3_open(dbfile.c_str(),&db);
  if( rc ) {
    fprintf(stderr,"Can't open database: %s\n",sqlite3_errmsg(db));
    return 0;
  }

  std::vector<dbEntry> entries(images.size());
  for(int q=0;q<images.size();q++){
    // Create SQL statement
    double k0 = images[q]["k"].asDouble();
    double g0 = images[q]["g"].asDouble();
    double s0 = images[q]["s"].asDouble();
    char buffer[250];
    sprintf(buffer,"SELECT id,k,g,s,(k-%f)*(k-%f)+(g-%f)*(g-%f) AS dkg,ABS(s-%f) AS ds,ss,res FROM gerlumph GROUP BY dkg,ds HAVING dkg < %f ORDER BY dkg ASC,ds ASC LIMIT 1;",k0,k0,g0,g0,s0,kg_sep*kg_sep);
    const char* sql = buffer;
    
    char count_buffer[250];
    sprintf(count_buffer,"SELECT COUNT(*) AS count,(k-%f)*(k-%f)+(g-%f)*(g-%f) AS dkg,ABS(s-%f) AS ds FROM gerlumph GROUP BY dkg,ds HAVING dkg < %f ORDER BY dkg ASC,ds ASC;",k0,k0,g0,g0,s0,kg_sep*kg_sep);
    const char* count_sql = count_buffer;


    // Execute SQL statements
    int row_count = 0;
    rc = sqlite3_exec(db,count_sql,count_callback,&row_count,&zErrMsg);
    if( rc != SQLITE_OK ) {
      fprintf(stderr,"SQL error: %s\n",zErrMsg);
      sqlite3_free(zErrMsg);
    }
    //fprintf(stdout,"%d\n",row_count);
    
    if( row_count == 0 && q != i_dtmax ){
      fprintf(stderr,"No matching map found for image %d within the given tolerance!\n",q);
      return 0;
    } else {
      dbEntry entry;
      const char* data = "Callback function called";
      rc = sqlite3_exec(db,sql,callback,&entry,&zErrMsg);
      if( rc != SQLITE_OK ) {
	fprintf(stderr,"SQL error: %s\n",zErrMsg);
	sqlite3_free(zErrMsg);
      }
      entries[q] = entry;
    }

  }
  sqlite3_close(db);


  // Write output JSON
  Json::Value json_db_entries;
  for(int i=0;i<entries.size();i++){
    Json::Value entry;
    entry["id"]   = entries[i].id;
    entry["k"]    = entries[i].k;
    entry["g"]    = entries[i].g;
    entry["s"]    = entries[i].s;
    entry["dkg"]  = std::sqrt(entries[i].dkg);
    entry["ds"]   = entries[i].ds;
    entry["width"] = entries[i].ss;
    entry["resolution"] = entries[i].res;
    json_db_entries.append(entry);
  }

  // JSON object needs to written like below because of double precision, e.g. 1.5 is written as 1.499999999999 otherwise
  std::ofstream file_maps(out_path+"output/gerlumph_maps.json");
  Json::StreamWriterBuilder wbuilder;
  wbuilder.settings_["precision"] = 6;
  std::unique_ptr<Json::StreamWriter> writer(wbuilder.newStreamWriter());
  writer->write(json_db_entries,&file_maps);
  file_maps.close();  
}
