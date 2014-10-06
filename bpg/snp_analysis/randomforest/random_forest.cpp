// mytest.cpp - displays a message

#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <string>
using namespace std;



class Node
{
public:
  Node(string ID, int depth, Node* parent);
  Node* GetParent(void);
  void CreateChildren(int maxdepth);
  void PrintNode(void);
  void GetChildren(int maxdepth, queue<Node*>* subtree);
private:
  string myID;
  int myDepth;
  Node* myParent; 
  Node* myLeftChild;
  Node* myRightChild;
};
Node::Node(string ID, int depth, Node* parent) {
  myID = ID;
  myDepth = depth;
  myParent = parent;
  myLeftChild = NULL;
  myRightChild = NULL;
}
Node* Node::GetParent(void) {
  return myParent;
}
void Node::CreateChildren(int maxdepth) {
  Node* left = new Node(this->myID+"1",this->myDepth+1,this);
  Node* right = new Node(this->myID+"2",this->myDepth+1,this);
  this->myLeftChild = left;  this->myRightChild = right;
  if (this->myDepth < maxdepth-1){ // check maxdepth not reached by children just created
    left->CreateChildren(maxdepth);
    right->CreateChildren(maxdepth);
  }
}
void Node::GetChildren(int maxdepth, queue<Node*>* subtree) {
  (*subtree).push(this);
  if (this->myLeftChild != NULL){ // intermediate nodes may not have children
    this->myLeftChild->GetChildren(maxdepth,subtree);
  }
  if (this->myRightChild != NULL){ // intermediate nodes may not have children
    this->myRightChild->GetChildren(maxdepth,subtree);
  }
}
void Node::PrintNode(void) {
  string parent, leftChild, rightChild;
  if (myParent != NULL){
    parent = myParent->myID;
  } else{
    parent = "None";
  }
  if (myLeftChild != NULL){
    leftChild = myLeftChild->myID;
  } else{
    leftChild = "None";
  }
  if (myRightChild != NULL){
    rightChild = myRightChild->myID;
  } else{
    rightChild = "None";
  }
  cout << myID << "\t\t" << myDepth << "\t\t" << parent << "\t\t" << leftChild << "\t\t" << rightChild << endl;
}


class Tree
{
public:
  Tree(string name, int maxdepth);
  string GetName(void);
  void GrowTree(void);
  void PrintTree(void);
private:
  string myName;
  int myMaxDepth;
  Node* myRootNode;
};
Tree::Tree(string name, int maxdepth) {
  myName = name;
  myMaxDepth = maxdepth;
  myRootNode = new Node("0",0,NULL);
}
string Tree::GetName(void) {
  return myName;
}
void Tree::GrowTree(void) {
  myRootNode->CreateChildren(myMaxDepth);
}
void Tree::PrintTree(void) {
  cout << endl;
  cout << "_________________________________________________________________________________________________" << endl;
  cout << "Tree: " << "\t" << myName << endl;
  cout << "Node Details: " << endl;
  cout << "ID" << "\t\t" << "Depth" << "\t\t" << "Parent" << "\t\t" << "LeftChild" << "\t" << "RightChild" << endl;
  cout << "----------" << "\t" << "----------" << "\t" << "----------" << "\t" << "----------" << "\t" << "----------" << endl;
  queue<Node*> subtree;
  myRootNode->GetChildren(myMaxDepth,&subtree);
  int subtreeSize = subtree.size();
  for (int i=0; i<subtreeSize; i++) {
    subtree.front()->PrintNode();
    subtree.pop();
  }
  cout << "_________________________________________________________________________________________________" << endl << endl;
}


class RandomForest
{
public:
  RandomForest(string name, int numtrees);
  string GetName(void);
  string GetTreeName(int treenum);
  int GetNumberTrees(void);
  void AddTree(Tree);
private:
  string myName;
  int myNumberOfTrees;
  vector<Tree> myTrees;
};

RandomForest::RandomForest(string name, int numtrees){
  myName = name;
  myNumberOfTrees = numtrees;
  vector<Tree*> myTrees;
}
string RandomForest::GetName(void){
  return myName;
}
string RandomForest::GetTreeName(int treenum){
  return myTrees[treenum].GetName();
}
int RandomForest::GetNumberTrees(void){
  return myNumberOfTrees;
}
void RandomForest::AddTree(Tree newtree){
  myTrees.push_back(newtree);
}


int main() {
  RandomForest forest("SomeForest",10);
  for (int i=0; i<forest.GetNumberTrees(); i++){
    Tree newTree("Some Tree",4);
    newTree.GrowTree();
    forest.AddTree(newTree);
    newTree.PrintTree();
  }
  //  for (int i=0; i<forest.GetNumberTrees(); i++){
  //cout << forest.GetTreeName(i) << endl;
  //}
  return 0;

}
