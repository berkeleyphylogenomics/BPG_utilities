import psycopg2

try:
    conn = psycopg2.connect("dbname='postgres' user='postgres' host='localhost' password='yaiyai'");
except:
    print "No DB"
    exit();

#DB won't update otherwise
conn.set_isolation_level(0)

cur = conn.cursor()

#cur.execute('INSERT INTO test_data(test1, test2) VALUES(66, 33)')

#cur.execute("SELECT * FROM test_data")

#rows = cur.fetchall()



