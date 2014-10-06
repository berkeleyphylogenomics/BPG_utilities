#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <mysql++/mysql++.h>
#include <mysql++/query.h>
#include <mysql++/result.h>
#include <mysql++/connection.h>

#define N 10
using namespace std;

typedef struct 
{
 long int id,parent_id,left_id,right_id;
 char scientific_name[256];
 char common_name[256];
}taxonomy;

typedef struct node 
{
 struct node *child[N];
 struct node *parent;
 taxonomy value;
}node;

const char *db="pfacts002",*server="phylodb",*user="rashmit",*pass="rhythm";
mysqlpp::Connection conn(false);
bool a = conn.connect(db,server,user,pass);
node *root,*start;

void insert_into_tree(long int id2,long int common_ancestor_id)
{
  long int parent_id2,left_id2,right_id2;
  char scientific_name2[256],common_name2[256];
  int i,j;

  if(id2 == common_ancestor_id)
  {
  	root = start;
	return;
  }
  mysqlpp::Query query1 = conn.query();	
  query1 << "SELECT id, parent_id, left_id, right_id, scientific_name, common_name  FROM ncbi_taxonomy WHERE ncbi_taxonomy.id = "<<id2;
  cout<<query1.str()<<endl; 
  mysqlpp::Result res4 = query1.store();
  cout << "Executed query " << endl;
  mysqlpp::Row row3 = res4.fetch_row();
  cout << "Fetched result " << endl;

  parent_id2 = row3["parent_id"];
  left_id2 = row3["left_id"];
  right_id2 = row3["right_id"];
  strcpy(scientific_name2,row3["scientific_name"]);
  strcpy(common_name2,row3["common_name"]);

  insert_into_tree(parent_id2,common_ancestor_id);
  
  for(i=0;i<=N;i++)
  {
	if(root->child[i] == NULL)
        {
	root->child[i] = (node *)malloc(sizeof(node));
	root->child[i]->value.id = id2;
	root->child[i]->value.parent_id = parent_id2;
	root->child[i]->value.left_id = left_id2;
	root->child[i]->value.right_id = right_id2;
	root->child[i]->parent = root;
	strcpy(root->child[i]->value.scientific_name,scientific_name2);
	strcpy(root->child[i]->value.common_name,common_name2);
	for(j=0;j<=N;j++)
		     root->child[i]->child[j] = NULL;
        root = root->child[i];
        goto ABC;
        }
	else if(root->child[i]->value.id == id2)
	{
        root = root->child[i];
	goto ABC;
	}
  }	
  ABC:
  return;
}

void create_NHX_tree(node *);
void print_node(node *);

ofstream outFile;

main()
{
 long int id1;
 long int min_left_id=20000000,max_right_id=0;
 long int left,right;
 long int common_ancestor_id;
 int i,j;
 int n;         // stores the no. of species as given in the input file

 ifstream inFile;
 inFile.open("taxonomyids.txt");
 if(!inFile){
 	cerr<<"Can't open input file "<<"taxonomyids.txt"<<endl;
	exit(1);
 }
 if(a)     //conn.connect(db,server,user,pass))
 {  
  {
  mysqlpp::NoExceptions ne(conn);
    i = 0;
    while(!inFile.eof())
    {
 	inFile>>id1;
	i++;
    }
    n = i-1;
    cout<<endl<<"No. of species : "<<n<<endl;
    inFile.close();
    inFile.open("taxonomyids.txt");
    for(i=0;i<n;i++)
    {
        inFile>>id1;
        cout<<endl<<id1<<" : ";   

        /* Find the minimum left id and maximum right id */
     try
     {
	mysqlpp::Query taxa_query = conn.query();
    	taxa_query<<"SELECT ncbi_taxonomy.left_id, ncbi_taxonomy.right_id FROM ncbi_taxonomy WHERE ncbi_taxonomy.id = "<<id1;
        mysqlpp::Result res = taxa_query.store();
	mysqlpp::Row row = res.fetch_row();
        
 		left = row["left_id"];
 		right = row["right_id"];
                
		cout<<"("<<left<<","<<right<<")" << endl;
 		if(left<min_left_id)
			min_left_id = left;
		if(right>max_right_id)
			max_right_id = right;
	//	i++;
     }
     catch(mysqlpp::Exception &e)
     {
        cout<<endl<<" BAD QUERY "<<endl;
	return false;
     }
    }
    //n = i;
    inFile.close();


   /* Code for finding the common ancestor */

   try {
   cout << "min_left_id: " << min_left_id << ", max_right_id: " << max_right_id;
   mysqlpp::Query new_query = conn.query();

   new_query << "SELECT id, parent_id, left_id, right_id, scientific_name, common_name FROM ncbi_taxonomy WHERE ncbi_taxonomy.left_id < "<<min_left_id<<" AND ncbi_taxonomy.right_id > "<<max_right_id<<" ORDER BY ncbi_taxonomy.left_id DESC LIMIT 1";
   cout<<new_query.str() << endl;
   mysqlpp::Result res1 = new_query.store();
   cout << "Executed query " << endl;
   mysqlpp::Row row1 = res1.fetch_row(); 
   cout << "Fetched result " << endl;
   

	common_ancestor_id = row1["id"];
   	root = (node *)malloc(sizeof(node));
   	root->value.id = common_ancestor_id;
	root->value.parent_id = row1["parent_id"];
	root->value.left_id = row1["left_id"];
	root->value.right_id = row1["right_id"];
	strcpy(root->value.scientific_name,row1["scientific_name"]);
	strcpy(root->value.common_name,row1["common_name"]);
	for(i=0;i<=N;i++)
		root->child[i] = NULL;
   
}
catch (const mysqlpp::BadQuery &er) {
	cerr << "Query error: " << er.what() << endl;
	exit(-1);
}
   
   /* Code for inserting each species into the tree */
  start = root;
  inFile.open("taxonomyids.txt",ios::in);
  for(i=0;i<n;i++)
  {
        inFile>>id1;
	insert_into_tree(id1,common_ancestor_id); 
   }
  inFile.close();
  create_NHX_tree(start);
  }
 }
 else
 {
 	cerr<<"Couldn't connect to the database"<<endl;
	exit(1);
 }
 conn.close();
}

void create_NHX_tree(node *root)
{
 outFile.open("species_tree.nhx",ios::out);
 if(!outFile)
 {
 	cerr<<"cannot open output file "<<endl;
	exit(1);
 }
 //outFile<<"(";
 print_node(root);
 //outFile<<")";
 outFile<<";";
 outFile.close();
}

/* This part is to print the species tree as required by NOTUNG, i.e. it removes those nodes in between which has only one child */

void print_node(node *root)
{
 int i;
 //outFile<<root->value.id;
 if(root->child[1] == NULL && root->child[0] != NULL)
 {
   print_node(root->child[0]);
   return;
 }
 if(root->child[0]!=NULL)
 {
	outFile<<"(";
   	i=0;
        while(root->child[i]!=NULL)
	{
	    print_node(root->child[i]);
	    if(root->child[i+1]!=NULL)
		outFile<<",";
 	    i=i+1;
	}
	outFile<<")";
 }
 outFile<<root->value.scientific_name;
}

/* This part is to print all the nodes in between, i.e. the species tree in whole containing all the taxons */ 

/*
	void print_node(node *root)
	{	
 	int i;
 	//outFile<<root->value.id;
 	if(root->child[0]!=NULL)
 	{
	outFile<<"(";
   	i=0;
	while(root->child[i]!=NULL)
	{
	    print_node(root->child[i]);
	    if(root->child[i+1]!=NULL)
		outFile<<",";
 	    i=i+1;
	}
	outFile<<")";
 	}
 	outFile<<root->value.scientific_name;
	}
*/
