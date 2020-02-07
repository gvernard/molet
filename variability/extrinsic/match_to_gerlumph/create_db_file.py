import sqlite3
 

conn = sqlite3.connect("gerlumph.db")
cur = conn.cursor()
cur.execute("CREATE TABLE gerlumph (id,k,g,s);") # use your column names here

# Using readline() 
file = open("/home/george/myCodes/molet/data/gerlumph_map_list/gerlumph_map_list.dat","r") 
  
while True: 
    line = file.readline() 
    if not line: 
        break
    dum = line.strip().split(" ")
    cur.execute("INSERT INTO gerlumph (id,k,g,s) VALUES (?,?,?,?);",dum)
file.close() 

conn.commit()
conn.close()
