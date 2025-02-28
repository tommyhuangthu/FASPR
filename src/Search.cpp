/*******************************************************************************************************************************
This file is a part of the protein side-chain packing software FASPR

Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/
#pragma warning(disable:4018)
#include "Search.h"


#define EPAIR_CUT      2.0
#define TREEWIDTH_CUT  5

void ShowGraph(Graph &graph){
  for(Graph::iterator it=graph.begin(); it!=graph.end(); ++it){
    cout <<it->first<<" =>";
    for(set<int>::iterator it2=it->second.begin();it2!=it->second.end();++it2){
      cout<<" "<<*it2;
    }
    cout<<endl;
  }
}


void Bag::ShowBag(){
  cout<<"(";
  for(set<int>::iterator it=left.begin();it!=left.end();++it){
    cout<<*it<<" ";
  }
  cout<<"|";
  for(set<int>::iterator it=right.begin();it!=right.end();++it){
    cout<<" "<<*it;
  }
  cout<<")";
}


void TreeDecomposition::Subgraph2TreeDecomposition(int index,Graph &graph){
  //sort the vertices in the cluster by the number of edges
  int ver=1;
  set <int> neibs;
  RankSize <set <int> > nbsize(ver, neibs);
  vector<RankSize< set<int> > > sortgraph;
  for(Graph::iterator it=graph.begin(); it!=graph.end(); ++it){
    ver=it->first;
    neibs=it->second;
    nbsize.idx=ver;
    nbsize.element=neibs;
    sortgraph.push_back(nbsize);
  }
  stable_sort(sortgraph.begin(), sortgraph.end(), less <RankSize<set<int> > >());
  cout<<endl<<"*subgraph "<<index<<" sorted (low-to-high degree):"<<endl;
  for(int i=0;i<sortgraph.size();++i){
    RankSize<set<int> > &st=sortgraph[i];
    cout <<st.idx<<" => ";
    for(set<int>::iterator it2=st.element.begin();it2!=st.element.end();++it2){
      cout<<" "<<*it2;
    }
    cout<<endl;
  }

  //construct the bags using a loop
  while (sortgraph.size()>0){
    //add to bag
    Bag newbag;
    RankSize< set<int> > &st=sortgraph[0];
    newbag.right.insert(st.idx);
    newbag.total.insert(st.idx);
    for(set<int>::iterator it2=st.element.begin();it2!=st.element.end();++it2){
      newbag.left.insert(*it2);
      newbag.total.insert(*it2);
    }
    bags.push_back(newbag);

    //1.remove edges from the sorted graph
    sortgraph.erase(sortgraph.begin());
    for(set<int>::iterator it3=newbag.left.begin();it3!=newbag.left.end();++it3){
      for(vector<RankSize<set<int> > >::iterator it4=sortgraph.begin();it4!=sortgraph.end();++it4){
        if(*it3==it4->idx){
          for(set<int>::iterator it5=it4->element.begin();it5!=it4->element.end();){
            if(*(newbag.right.begin())==*it5){
              it4->element.erase(it5++);
              if(it4->element.empty()){
                break;
              }
            }
            else{
              it5++;
            }
          }
          if(it4->element.empty()){
            sortgraph.erase(it4++);
            break;
          }
        }
      }
    }

    //2.add edges for any two residues in the left bag
    if(newbag.left.size()>1){
      for(set<int>::iterator it3=newbag.left.begin();it3!=newbag.left.end();++it3){
        bool exist=false;
        for(int it5=0;it5<sortgraph.size();++it5){
          if(sortgraph[it5].idx==*it3){
            exist=true;
            break;
          }
        }
        if(exist==false){
          set <int> neibs;
          RankSize< set<int> > newedge(*it3,neibs);
          sortgraph.push_back(newedge);
        }
      }
      for(set<int>::iterator it3=newbag.left.begin();it3!=newbag.left.end();++it3){
        for(set<int>::iterator it4=it3;it4!=newbag.left.end();++it4){
          if(it4==it3) continue;
          for(int it5=0;it5<sortgraph.size();++it5){
            if(sortgraph[it5].idx==*it3){
              sortgraph[it5].element.insert(*it4);
              break;
            }
          }

          for(int it5=0;it5<sortgraph.size();++it5){
            if(sortgraph[it5].idx==*it4){
              sortgraph[it5].element.insert(*it3);
              break;
            }
          }
        }
      }
    }
    
    //3.resort the subgraph
    stable_sort(sortgraph.begin(), sortgraph.end(), less <RankSize<set<int> > >());
  }

  //show all of the bags on the tree
  if(true){
    cout<<"@list bags (nodes) of the tree decomposition:"<<endl;
    for(int i=0; i<bags.size(); ++i){
      bags[i].ShowBag();
      cout<<endl;
    }
  }

  //connect the bags into a tree
  set<int> vtsOnTree;
  int counter=0;
  while(bags.size()>0){
    Bag bg=bags[bags.size()-1];
    bags.pop_back();
    if(counter==0){
      //for the root bag, set the parent Bag Index as -1
      bg.parentBagIdx=-1;
      bg.type=Root;
      connBags.push_back(bg);
      //add vertices from bag to vtsOnTree;
      set_union(vtsOnTree.begin(),vtsOnTree.end(),
        bg.total.begin(),bg.total.end(),
        inserter(vtsOnTree,vtsOnTree.begin())
        );
    }
    else{
      set<int> intersection;
      set_intersection(vtsOnTree.begin(),vtsOnTree.end(),
        bg.total.begin(),bg.total.end(),
        inserter(intersection,intersection.begin())
        );
      for(int i=0; i<connBags.size();++i){
        Bag nbg=connBags[i];
        set<int> intersection2;
        set_intersection(nbg.total.begin(),nbg.total.end(),
          intersection.begin(),intersection.end(),
          inserter(intersection2,intersection2.begin())
          );
        if(intersection2==intersection){
          //bag nbg is the parent of the current bag
          //modify the attributes of parent bag
          bg.parentBagIdx=i;
          if(connBags[i].type!=Root){
            connBags[i].type=Inner;
          }
          connBags[i].childBagIdx.insert(connBags.size());
          connBags.push_back(bg);
          //add vertices from bag to vtsOnTree;
          set_union(vtsOnTree.begin(),vtsOnTree.end(),
            bg.total.begin(),bg.total.end(),
            inserter(vtsOnTree,vtsOnTree.begin())
            );
          break;
        }
      }
    }
    counter++;
  }
  //set the number of children for each bag
  for(int i=0; i<connBags.size(); ++i){
    connBags[i].childCounter=connBags[i].childBagIdx.size();
  }

  //show the connected bags, parent and children
  if(false){
    cout<<"@list connected bags of the tree decomposition:"<<endl;
    for(int i=0; i<connBags.size(); ++i){
      connBags[i].ShowBag();
      if(connBags[i].type==Root){
        cout<<":"<<endl<<"index = "<<setw(4)<<i<<", type =  root, parent = none";
        cout<<", children =";
        for(set<int>::iterator it1=connBags[i].childBagIdx.begin(); it1!=connBags[i].childBagIdx.end();++it1){
          cout<<" "<<*it1;
        }
        cout<<endl;
      }
      else if(connBags[i].type==Inner){
        cout<<":"<<endl<<"index = "<<setw(4)<<i<<", type = inner, parent = "<<setw(4)<<connBags[i].parentBagIdx;
        cout<<", children =";
        for(set<int>::iterator it1=connBags[i].childBagIdx.begin(); it1!=connBags[i].childBagIdx.end();++it1){
          cout<<" "<<*it1;
        }
        cout<<endl;
      }
      else{
        cout<<":"<<endl<<"index = "<<setw(4)<<i<<", type =  leaf, parent = "<<setw(4)<<connBags[i].parentBagIdx;
        cout<<", children = none"<<endl;
      }
    }
  }

}

void TreeDecomposition::MergeBags(int depth)
{
  if(depth<connBags.size()){
    for(int i=0;i<connBags.size();i++){
      if(connBags[i].type==Leaf){
        //merge the leaf to its parent
        cout<<"Merge leaf "<<i<<" with its parent "<<connBags[i].parentBagIdx<<": ";
        set_union(connBags[i].total.begin(),connBags[i].total.end(),
          connBags[connBags[i].parentBagIdx].total.begin(),connBags[connBags[i].parentBagIdx].total.end(),
          inserter(connBags[connBags[i].parentBagIdx].total,connBags[connBags[i].parentBagIdx].total.begin())
          );
        for(set<int>::iterator it=connBags[connBags[i].parentBagIdx].total.begin();it!=connBags[connBags[i].parentBagIdx].total.end();++it){
          cout<<*it<<" ";
        }
        cout<<endl;
        connBags[i].type=None;
        connBags[connBags[i].parentBagIdx].childBagIdx.erase(i);
        if(connBags[connBags[i].parentBagIdx].childBagIdx.size()==0 && connBags[connBags[i].parentBagIdx].type != Root){
          connBags[connBags[i].parentBagIdx].type=Leaf;
        }
        depth++;
      }
    }
    MergeBags(depth);
  }
}

int TreeDecomposition::CheckTreewidth()
{
  int width=0;
  for(int i=0;i<connBags.size();i++){
    if(width<connBags[i].total.size()){
      width=connBags[i].total.size();
    }
  }
  width++;
  return width;
}


Solution::~Solution()
{
  unfixres.clear();
}


bool Solution::DEESearch(IV1 &pos){
  //DEEGoldstein
  IV1 fixres;
  int ndeadDEE=1,iterDEE=1;
  unfixres.clear();
  while(ndeadDEE!=0){
    ndeadDEE=DEEGoldstein(pos);
    cout<<"iter "<<iterDEE<<" DEEgoldstein eliminates "<<ndeadDEE<<" rotamers"<<endl;
    iterDEE++;
  };
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    int n=0;
    int rot;
    for(int j=0;j<nrots[ip];j++){
      if(eTableSelf[ip][j]<999.){
        rot=j;
        n++;
      }
    }
    if(n==1){
      nrots[ip]=1;
      bestrot[ip]=rot;
      fixres.push_back(ip);
    }
    else if(n>1){
      unfixres.push_back(ip);
    }
  }
  cout<<"#residues fixed after DEE-Goldstein: "<<fixres.size()<<endl;
  cout<<"#residues unfixed after DEE-Goldstein: "<<unfixres.size()<<endl;
  if(unfixres.size()==0) return false;
  for(int i=0;i<unfixres.size();i++){
    int ipos=unfixres[i];
    for(int j=0;j<fixres.size();j++){
      int jpos=fixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      int rot=bestrot[jpos];
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]<999.){
          eTableSelf[ipos][k] += eTablePair[ipos][jpos][k][rot];
        }
      }
    }
  }

  //DEEsplit
  pos.clear();
  fixres.clear();
  unfixres.clear();
  ndeadDEE=1,iterDEE=1;
  for(int i=0;i<nres;i++){
    if(nrots[i]>1) pos.push_back(i);
  }
  while(ndeadDEE!=0){
    ndeadDEE=DEEsplit(pos);
    cout<<"iter "<<iterDEE<<" DEEsplit eliminates "<<ndeadDEE<<" rotamers"<<endl;
    iterDEE++;
  };
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    int n=0;
    int rot;
    for(int j=0;j<nrots[ip];j++){
      if(eTableSelf[ip][j]<999.){
        rot=j;
        n++;
      }
    }
    if(n==1){
      nrots[ip]=1;
      bestrot[ip]=rot;
      fixres.push_back(ip);
    }
    else if(n>1){
      unfixres.push_back(ip);
    }
  }
  cout<<"#residues fixed after DEE-split: "<<fixres.size()<<endl;
  cout<<"#residues unfixed after DEE-split: "<<unfixres.size()<<endl;
  if(unfixres.size()==0) return false;
  for(int i=0;i<unfixres.size();i++){
    int ipos=unfixres[i];
    for(int j=0;j<fixres.size();j++){
      int jpos=fixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      int rot=bestrot[jpos];
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]<999.){
          eTableSelf[ipos][k] += eTablePair[ipos][jpos][k][rot];
        }
      }
    }
  }

  //DEEGoldstein
  pos.clear();
  fixres.clear();
  unfixres.clear();
  ndeadDEE=1,iterDEE=1;
  for(int i=0;i<nres;i++){
    if(nrots[i]>1) pos.push_back(i);
  }
  while(ndeadDEE!=0){
    ndeadDEE=DEEGoldstein(pos);
    cout<<"iter "<<iterDEE<<" DEEgoldstein eliminates "<<ndeadDEE<<" rotamers"<<endl;
    iterDEE++;
  };
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    int n=0;
    int rot;
    for(int j=0;j<nrots[ip];j++){
      if(eTableSelf[ip][j]<999.){
        rot=j;
        n++;
      }
    }
    if(n==1){
      nrots[ip]=1;
      bestrot[ip]=rot;
      fixres.push_back(ip);
    }
    else if(n>1){
      unfixres.push_back(ip);
    }
  }
  cout<<"#residues fixed after DEE-Goldstein: "<<fixres.size()<<endl;
  cout<<"#residues unfixed after DEE-Goldstein: "<<unfixres.size()<<endl;
  if(unfixres.size()==0) return false;
  for(int i=0;i<unfixres.size();i++){
    int ipos=unfixres[i];
    for(int j=0;j<fixres.size();j++){
      int jpos=fixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      int rot=bestrot[jpos];
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]<999.){
          eTableSelf[ipos][k] += eTablePair[ipos][jpos][k][rot];
        }
      }
    }
  }

  //DEEsplit
  pos.clear();
  fixres.clear();
  unfixres.clear();
  ndeadDEE=1,iterDEE=1;
  for(int i=0;i<nres;i++){
    if(nrots[i]>1) pos.push_back(i);
  }
  while(ndeadDEE!=0){
    ndeadDEE=DEEsplit(pos);
    cout<<"iter "<<iterDEE<<" DEEsplit eliminates "<<ndeadDEE<<" rotamers"<<endl;
    iterDEE++;
  };
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    int n=0;
    int rot;
    for(int j=0;j<nrots[ip];j++){
      if(eTableSelf[ip][j]<999.){
        rot=j;
        n++;
      }
    }
    if(n==1){
      nrots[ip]=1;
      bestrot[ip]=rot;
      fixres.push_back(ip);
    }
    else if(n>1){
      unfixres.push_back(ip);
    }
  }
  cout<<"#residues fixed after DEE-split: "<<fixres.size()<<endl;
  cout<<"#residues unfixed after DEE-split: "<<unfixres.size()<<endl;
  if(unfixres.size()==0) return false;
  for(int i=0;i<unfixres.size();i++){
    int ipos=unfixres[i];
    for(int j=0;j<fixres.size();j++){
      int jpos=fixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      int rot=bestrot[jpos];
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]<999.){
          eTableSelf[ipos][k] += eTablePair[ipos][jpos][k][rot];
        }
      }
    }
  }

  return true;
}

int Solution::DEEGoldstein(IV1& pos)
{
  int elimination=0;
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    for(int s=0;s<nrots[ip];s++){
      if(eTableSelf[ip][s]>999.) continue;
      for(int r=0;r<nrots[ip];r++){
        if(r==s) continue;
        else if(eTableSelf[ip][r]>999.) continue;
        float ex=eTableSelf[ip][s]-eTableSelf[ip][r];
        for(int j=0;j<pos.size();j++){
          if(j==i) continue;
          int jp=pos[j];
          if(eTablePair[ip][jp]==NULL) continue;
          float ey=1e8;
          for(int t=0;t<nrots[jp];t++){
            if(eTableSelf[jp][t]>999.) continue;
            float em=eTablePair[ip][jp][s][t]-eTablePair[ip][jp][r][t];
            if(em<ey) ey=em;
          }
          ex+=ey;
        }
        if(ex>DEE_THRESHOLD){
          eTableSelf[ip][s]=1000.;
          elimination++;
          break;
        }
      }
    }
  }
  return elimination;
}

int Solution::DEEsplit(IV1& pos)
{
  int elimination=0;
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    for(int s=0;s<nrots[ip];s++){
      if(eTableSelf[ip][s]>999.) continue;
      FV1 storeYj(pos.size(),0.);
      FV2 storeYjr(nrots[ip],storeYj);
      for(int r=0;r<nrots[ip];r++){
        if(r==s) continue;
        else if(eTableSelf[ip][r]>999.) continue;
        for(int j=0;j<pos.size();j++){
          if(j==i) continue;
          int jp=pos[j];
          if(eTablePair[ip][jp]==NULL) continue;
          float ey=1e8;
          for(int t=0;t<nrots[jp];t++){
            if(eTableSelf[jp][t]>999.) continue;
            float em=eTablePair[ip][jp][s][t]-eTablePair[ip][jp][r][t];
            if(em<ey) ey=em;
          }
          storeYjr[r][j]=ey;
        }
      }

      IV1 elim;
      for(int k=0;k<pos.size();k++){
        if(k==i) continue;
        int kp=pos[k];
        if(eTablePair[ip][kp]==NULL) continue;
        for(int v=0;v<nrots[kp];v++){
          elim.push_back(0);
        }
        for(int r=0;r<nrots[ip];r++){
          if(r==s) continue;
          else if(eTableSelf[ip][r]>999.) continue;

          float ex=eTableSelf[ip][s]-eTableSelf[ip][r];
          for(int j=0;j<pos.size();j++){
            if(j==i || j==k) continue;
            int jp=pos[j];
            if(eTablePair[ip][jp]==NULL) continue;
            ex+=storeYjr[r][j];
          }
          for(int v=0;v<nrots[kp];v++){
            if(eTableSelf[kp][v]>999.) continue;
            if(ex+eTablePair[ip][kp][s][v]-eTablePair[ip][kp][r][v]>DEE_THRESHOLD){
              elim[v]=1;
            }
          }
        }
        bool allelim=true;
        for(int v=0;v<nrots[kp];v++){
          if(elim[v]==0){
            allelim=false;
            break;
          }
        }
        if(allelim==true){
          eTableSelf[ip][s]=1000.;
          elimination++;
          goto FLAG_SPLIT;
        }
      }
FLAG_SPLIT:;
    }
  }
  return elimination;
}

void Solution::Pick(int site,int rot)
{
  for(int i=5;i<stru[site].atTypes.size();i++){
    stru[site].xyz.push_back(sc[site][rot][i-5]);
  }
}

void Solution::ConstructAdjMatrix(int nunfix,IV2 &adjMatrix)
{
  cout<<endl<<"construct adjacent matrix ..."<<endl;
  cout<<"remove edges (residue-residue interactions) with small energy values ... ";
  int n_removed_edges=0;
  for(int i=0;i<nunfix-1;i++){
    int ipos=unfixres[i];
    for(int j=i+1;j<nunfix;j++){
      int jpos=unfixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]>999.)continue;
        for(int l=0;l<nrots[jpos];l++){
          if(eTableSelf[jpos][l]>999.)continue;
          if(eTablePair[ipos][jpos][k][l]>EPAIR_CUT || eTablePair[ipos][jpos][k][l]<-1.*EPAIR_CUT){
            adjMatrix[i][j]=1;
            adjMatrix[j][i]=1;
            goto FLAG1;
          }
        }
      }
      eTablePair[ipos][jpos]=NULL;
      //cout<<"residue "<<ipos<<" and "<<jpos<<" does not interact"<<endl;
      n_removed_edges++;
FLAG1: continue;
    }
  }
  cout<<"#edges removed: "<<n_removed_edges<<endl;

  cout<<"remove the residues that have no interaction edge ... ";
  int n_res_noedge=0;
  for(int i=0;i<nunfix;i++){
    bool allzeros=true;
    for(int j=0;j<nunfix;j++){
      if(adjMatrix[i][j]==1){
        allzeros=false;
        break;
      }
    }
    if(allzeros==true){
      float emin=1e8;
      int rot;
      for(int r=0;r<nrots[unfixres[i]];r++){
        if(eTableSelf[unfixres[i]][r]<emin){
          emin=eTableSelf[unfixres[i]][r];
          rot=r;
        }
      }
      nrots[unfixres[i]]=1;
      bestrot[unfixres[i]]=rot;
      n_res_noedge++;
    }
  }
  cout<<"#residues without edges removed: "<<n_res_noedge<<endl;
}




void Solution::ConstructSubgraphs(int nunfix,IV1 &visited,IV2 &adjMatrix,IV2 &flagMatrix)
{
  visited.assign(nunfix,0);
  for(int i=0;i<nunfix;i++){
    if(!visited[i]){
      Graph newsg;
      stack<int> vertices;
      FindSubgraphByDFS(newsg,i,visited,adjMatrix,flagMatrix,vertices);
      if(!newsg.empty()){
        graphs.push_back(newsg);
      }

      if(newsg.size()==1){
        int site=unfixres[i];
        int rot=0;
        float eMin=1000.;
        for(int j=0;j<nrots[site];j++){
          if(eTableSelf[site][j]>999.)continue;
          if(eTableSelf[site][j]<eMin){
            eMin=eTableSelf[site][j];
            rot=j;
          }
        }
        bestrot[site]=rot;
        nrots[site]=1;
      }
    }
  }

  for(vector<int>::iterator it=unfixres.begin();it!=unfixres.end();){
    if(nrots[*it]==1){
      it=unfixres.erase(it);
    }
    else{
      ++it;
    }
  }


}


void Solution::FindSubgraphByDFS(Graph &graph,int u,IV1 &visited,IV2 &adjMatrix,IV2 &flagMatrix,stack<int> &vertices)
{
  visited[u]=1;
  vertices.push(u);
  while(!vertices.empty()){
    bool hasEdge=false;
    int w=-1;
    for(w=0;w<adjMatrix[u].size();++w){
      if(adjMatrix[u][w]==1 && flagMatrix[u][w]==0){
        hasEdge=true;
        break;
      }
    }

    if(hasEdge==false){
      u=vertices.top();
      vertices.pop();
    }
    else{
      visited[w]=1;
      vertices.push(w);
      flagMatrix[u][w]=1;
      flagMatrix[w][u]=1;
      graph[unfixres[u]].insert(unfixres[w]);
      graph[unfixres[w]].insert(unfixres[u]);
    }
  }
}


void Solution::ShowGraphs(){
  for(int i1=0; i1<graphs.size(); i1++){
    cout<<"*subgraph "<<i1<<":"<<endl;
    Graph &graph=graphs[i1];
    ShowGraph(graph);
  }
}



void Solution::GetLeftBagRotamerCombination(Bag &leafbag,int depth,IV1 &Rtmp,IV2 &Rlcom){
  if(depth<leafbag.lsites.size()){
    for(int i=0;i<leafbag.lrots[depth].size();i++){
      int roti=leafbag.lrots[depth][i];
      Rtmp.push_back(roti);
      GetLeftBagRotamerCombination(leafbag,depth+1,Rtmp,Rlcom);
      Rtmp.pop_back();
    }
  }
  else{
    Rlcom.push_back(Rtmp);
  }
}

void Solution::CalcLeftBagRotamerCombinationEnergy(Bag &rootbag,int depth,float &Etmp,IV1 &Rtmp,FV1 &Elcom,IV2 &Rlcom){
  if(depth<rootbag.lsites.size()){
    int site=rootbag.lsites[depth];
    for(int i=0;i<rootbag.lrots[depth].size();i++){
      int roti=rootbag.lrots[depth][i];
      float Eold=Etmp;
      float Enew=eTableSelf[site][roti];
      for(int k=0;k<depth;k++){
        if(eTablePair[rootbag.lsites[k]][site]==NULL) continue;
        Enew+=eTablePair[rootbag.lsites[k]][site][Rlcom[Rlcom.size()-1][k]][roti];
      }
      Etmp+=Enew;
      Rtmp.push_back(roti);
      CalcLeftBagRotamerCombinationEnergy(rootbag,depth+1,Etmp,Rtmp,Elcom,Rlcom);
      Rtmp.pop_back();
      Etmp=Eold;
    }
  }
  else{
    Rlcom.push_back(Rtmp);
    Elcom.push_back(Etmp);
  }
}



void Solution::CalcRightBagRotamerCombinationEnergy(Bag &leafbag,int depth,float &Etmp,IV1 &Rtmp,FV1 &Ercom,IV2 &Rrcom){
  if(depth<leafbag.rsites.size()){
    int site=leafbag.rsites[depth];
    for(int i=0;i<leafbag.rrots[depth].size();i++){
      int roti=leafbag.rrots[depth][i];
      float Eold=Etmp;
      float Enew=eTableSelf[site][roti];
      for(int k=0;k<depth;k++){
        if(eTablePair[leafbag.rsites[k]][site]==NULL) continue;
        Enew+=eTablePair[leafbag.rsites[k]][site][Rrcom[Rrcom.size()-1][k]][roti];
      }
      Etmp+=Enew;
      Rtmp.push_back(roti);
      CalcRightBagRotamerCombinationEnergy(leafbag,depth+1,Etmp,Rtmp,Ercom,Rrcom);
      Rtmp.pop_back();
      Etmp=Eold;
    }
  }
  else{
    Rrcom.push_back(Rtmp);
    Ercom.push_back(Etmp);
  }
}


void Solution::BagDeploySites(Bag &leafbag)
{
  if(leafbag.deployFlag==false){
    leafbag.deployFlag=true;
    for(set<int>::iterator it=leafbag.left.begin();it!=leafbag.left.end();it++){
      leafbag.lsites.push_back(*it);
      leafbag.tsites.push_back(*it);
    }
    for(set<int>::iterator it=leafbag.right.begin();it!=leafbag.right.end();it++){
      leafbag.rsites.push_back(*it);
      leafbag.tsites.push_back(*it);
    }

    IV1 rots;
    for(int idx2=0;idx2<leafbag.lsites.size();idx2++){
      int site=leafbag.lsites[idx2];
      for(int j=0;j<nrots[site];j++){
        if(eTableSelf[site][j]>999.) continue;
        rots.push_back(j);
      }
      leafbag.lrots.push_back(rots);
      rots.clear();
    }
    for(int idx2=0;idx2<leafbag.rsites.size();idx2++){
      int site=leafbag.rsites[idx2];
      for(int j=0;j<nrots[site];j++){
        if(eTableSelf[site][j]>999.) continue;
        rots.push_back(j);
      }
      leafbag.rrots.push_back(rots);
      rots.clear();
    }

    int lcount=1;
    int counter=0;
    while(counter<leafbag.lrots.size()){
      lcount *=leafbag.lrots[counter].size();
      counter++;
    }

    int rcount=1;
    counter=0;
    while(counter<leafbag.rrots.size()){
      rcount *= leafbag.rrots[counter].size();
      counter++;
    }
    for(int i=0;i<lcount;i++){
      IV1 Rtmp(rcount,-1);
      FV1 Etmp(rcount,0.);
      leafbag.Etcom.push_back(Etmp);
      leafbag.Rtcom.push_back(Rtmp);
    }
  }
}



void Solution::LeafBagCalcEnergy(Bag &leafbag,IV2 &Rlcom)
{
  FV1 Ercom;
  float Ertmp=0.;
  IV2 Rrcom;
  IV1 Rrtmp, Rltmp;
  CalcRightBagRotamerCombinationEnergy(leafbag,0,Ertmp,Rrtmp,Ercom,Rrcom);
  GetLeftBagRotamerCombination(leafbag,0,Rltmp,Rlcom);
  for(int idx2=0;idx2<Rlcom.size();idx2++){
    float emin=1e8;
    IV1 Rrmin;
    for(int idx3=0;idx3<Rrcom.size();idx3++){
      float eval=0.;
      for(int idx4=0;idx4<Rlcom[idx2].size();idx4++){
        for(int idx5=0;idx5<Rrcom[idx3].size();idx5++){
          if(eTablePair[leafbag.lsites[idx4]][leafbag.rsites[idx5]] != NULL){
            eval+=eTablePair[leafbag.lsites[idx4]][leafbag.rsites[idx5]][Rlcom[idx2][idx4]][Rrcom[idx3][idx5]];
          }
        }
      }
      eval+=Ercom[idx3];
      if(leafbag.Etcom.size() != 0){
        //add the pre-stored energy
        eval+=leafbag.Etcom[idx2][idx3];
      }
      if(eval<emin){
        Rrmin=Rrcom[idx3];
        emin=eval;
      }
    }
    leafbag.Etcom[idx2][0]=emin;

    //record the best rotamer combinations
    IV1 Rttmp;
    for(int idx4=0;idx4<Rlcom[idx2].size();idx4++){
      Rttmp.push_back(Rlcom[idx2][idx4]);
    }
    for(int idx5=0;idx5<Rrmin.size();idx5++){
      Rttmp.push_back(Rrmin[idx5]);
    }
    leafbag.Rtcom[idx2]=Rttmp;
    Rttmp.clear();
  }

  Ercom.clear();
  Rrcom.clear();
  Rrtmp.clear();
  Rltmp.clear();
}


void Solution::CombineChildIntoParentBag(Bag &leafbag,Bag &parbag,IV2 &Rclcom)
{
  FV1 Ercom;
  float Ertmp=0.;
  IV2 Rrcom, Rlcom;
  IV1 Rrtmp, Rltmp;
  CalcRightBagRotamerCombinationEnergy(parbag,0,Ertmp,Rrtmp,Ercom,Rrcom);
  GetLeftBagRotamerCombination(parbag,0,Rltmp,Rlcom);
  SubsetCheck(leafbag.lsites,parbag.tsites,leafbag.indices);
  for(int idx2=0;idx2<Rlcom.size();idx2++){
    for(int idx3=0;idx3<Rrcom.size();idx3++){
      IV1 ppartrot,ppartsite;
      for(int j=0;j<leafbag.indices.size();j++){
        if(leafbag.indices[j]>=parbag.lsites.size()){
          ppartrot.push_back(Rrcom[idx3][leafbag.indices[j]-parbag.lsites.size()]);
          ppartsite.push_back(parbag.rsites[leafbag.indices[j]-parbag.lsites.size()]);
        }
        else{
          ppartrot.push_back(Rlcom[idx2][leafbag.indices[j]]);
          ppartsite.push_back(parbag.lsites[leafbag.indices[j]]);
        }
      }

      for(int j=0;j<Rclcom.size();j++){
        if(ppartsite==leafbag.lsites && ppartrot==Rclcom[j]){
          parbag.Etcom[idx2][idx3]+=leafbag.Etcom[j][0];
        }
      }
    }
  }
}

void Solution::SubsetCheck(IV1 &leaflsites,IV1 &partsites, IV1 &indices)
{
  indices.clear();
  for(int i=0;i<leaflsites.size();i++){
    indices.push_back(-1);
    for(int j=0;j<partsites.size();j++){
      if(leaflsites[i]==partsites[j]){
        indices[i]=j;
      }
    }
  }
}



void Solution::RootBagFindGMEC(Bag &rootbag)
{
  FV1 Ercom,Elcom;
  float Ertmp=0., Eltmp=0.;
  IV2 Rrcom, Rlcom;
  IV1 Rrtmp, Rltmp;
  CalcRightBagRotamerCombinationEnergy(rootbag,0,Ertmp,Rrtmp,Ercom,Rrcom);
  CalcLeftBagRotamerCombinationEnergy(rootbag,0,Eltmp,Rltmp,Elcom,Rlcom);

  rootbag.EGMEC=1e8;
  for(int idx2=0;idx2<Rlcom.size();idx2++){
    for(int idx3=0;idx3<Rrcom.size();idx3++){

      float energy=Elcom[idx2];
      for(int idx4=0;idx4<Rlcom[idx2].size();idx4++){
        for(int idx5=0;idx5<Rrcom[idx3].size();idx5++){
          if(eTablePair[rootbag.lsites[idx4]][rootbag.rsites[idx5]] != NULL){
            energy+=eTablePair[rootbag.lsites[idx4]][rootbag.rsites[idx5]][Rlcom[idx2][idx4]][Rrcom[idx3][idx5]];
          }
        }
      }
      energy+=Ercom[idx3];
      if(rootbag.Etcom.size() != 0){
        energy+=rootbag.Etcom[idx2][idx3];
      }

      if(energy<rootbag.EGMEC){
        rootbag.EGMEC=energy;
        rootbag.RLGMEC.assign(Rlcom[idx2].begin(),Rlcom[idx2].end());
        rootbag.RRGMEC.assign(Rrcom[idx3].begin(),Rrcom[idx3].end());
      }

    }
  }

  //set optimal rotamer index
  for(int i=0;i<rootbag.lsites.size();i++){
    nrots[rootbag.lsites[i]]=1;
    bestrot[rootbag.lsites[i]]=rootbag.RLGMEC[i];
  }
  for(int i=0;i<rootbag.rsites.size();i++){
    nrots[rootbag.rsites[i]]=1;
    bestrot[rootbag.rsites[i]]=rootbag.RRGMEC[i];
  }

}

void Solution::TreeDecompositionBottomToTopCalcEnergy()
{
  if(tree.connBags.size()==1){
    BagDeploySites(tree.connBags[0]);
  }
  else{
    while(true){
      int leafCount=0;
      for(int idx1=0;idx1<tree.connBags.size();idx1++){
        if(tree.connBags[idx1].type==Leaf){
          IV2 Rclcom;
          BagDeploySites(tree.connBags[idx1]);
          LeafBagCalcEnergy(tree.connBags[idx1],Rclcom);

          BagDeploySites(tree.connBags[tree.connBags[idx1].parentBagIdx]);
          CombineChildIntoParentBag(tree.connBags[idx1],tree.connBags[tree.connBags[idx1].parentBagIdx],Rclcom);

          tree.connBags[idx1].type=None;
          tree.connBags[tree.connBags[idx1].parentBagIdx].childCounter--;
          if(tree.connBags[tree.connBags[idx1].parentBagIdx].childCounter==0 && 
            tree.connBags[tree.connBags[idx1].parentBagIdx].type != Root){
              tree.connBags[tree.connBags[idx1].parentBagIdx].type=Leaf;
          }
          leafCount++;
        }
      }
      if(leafCount==0){
        break;
      }
    }
  }

}


void Solution::TreeDecompositionTopToBottomAssignRotamer(Bag &parbag,Bag &childbag)
{
  bool find=false;
  IV1 ppartrot,ppartsite;
  for(int j=0;j<childbag.indices.size();j++){
    if(childbag.indices[j]>=parbag.lsites.size()){
      ppartrot.push_back(parbag.RRGMEC[childbag.indices[j]-parbag.lsites.size()]);
      ppartsite.push_back(parbag.rsites[childbag.indices[j]-parbag.lsites.size()]);
    }
    else{
      ppartrot.push_back(parbag.RLGMEC[childbag.indices[j]]);
      ppartsite.push_back(parbag.lsites[childbag.indices[j]]);
    }
  }

  for(int j=0;j<childbag.Rtcom.size();j++){
    IV1 Rcltmp,Rcrtmp;
    Rcltmp.assign(childbag.Rtcom[j].begin(),childbag.Rtcom[j].begin()+childbag.lsites.size());
    Rcrtmp.assign(childbag.Rtcom[j].begin()+childbag.lsites.size(),childbag.Rtcom[j].end());
    if(ppartsite==childbag.lsites && ppartrot==Rcltmp){
      childbag.RLGMEC=Rcltmp;
      childbag.RRGMEC=Rcrtmp;
      for(int p=0;p<childbag.rsites.size();p++){
        nrots[childbag.rsites[p]]=1;
        bestrot[childbag.rsites[p]]=childbag.RRGMEC[p];
      }
      find=true;
      break;
    }
  }

  for(set<int>::iterator it=childbag.childBagIdx.begin();it!=childbag.childBagIdx.end();it++){
    TreeDecompositionTopToBottomAssignRotamer(childbag,tree.connBags[*it]);
  }
}


void Solution::TreeDecompositionRelease()
{
  tree.bags.clear();
  tree.connBags.clear();
}


void Solution::GraphEdgeDecomposition(IV2 &adjMatrix,float threshold)
{
  int nedgeremoved=0;
  for(int i=0;i<adjMatrix.size()-1;i++){
    int k=unfixres[i];
    for(int j=i+1;j<adjMatrix.size();j++){
      int l=unfixres[j];
      if(adjMatrix[i][j]==1 && eTablePair[k][l]!=NULL){
        int counterm=0, countern=0;
        for(int m=0;m<nrots[k];m++){
          if(eTableSelf[k][m]>999.0) continue;
          counterm++;
        }
        for(int n=0;n<nrots[l];n++){
          if(eTableSelf[l][n]>999.0) continue;
          countern++;
        }

        float abar=0.;
        for(int m=0;m<nrots[k];m++){
          if(eTableSelf[k][m]>999.0) continue;
          for(int n=0;n<nrots[l];n++){
            if(eTableSelf[l][n]>999.0) continue;
            abar+=eTablePair[k][l][m][n];
          }
        }
        abar/=(2.*counterm*countern);

        FV1 ak,bl;
        ak.assign(nrots[k],1000.);
        bl.assign(nrots[l],1000.);
        for(int m=0;m<nrots[k];m++){
          float temp=0.;
          if(eTableSelf[k][m]>999.0) continue;
          for(int n=0;n<nrots[l];n++){
            if(eTableSelf[l][n]>999.0) continue;
            temp+=eTablePair[k][l][m][n];
          }
          temp/=countern;
          ak[m]=temp-abar;
        }

        for(int n=0;n<nrots[l];n++){
          float temp=0.;
          if(eTableSelf[l][n]>999.0) continue;
          for(int m=0;m<nrots[k];m++){
            if(eTableSelf[k][m]>999.0) continue;
            temp+=eTablePair[k][l][m][n];
          }
          temp/=counterm;
          bl[n]=temp-abar;
        }

        //estimate max deviation
        float maxdev=-1e8;
        for(int m=0;m<nrots[k];m++){
          if(eTableSelf[k][m]>999.0) continue;
          for(int n=0;n<nrots[l];n++){
            if(eTableSelf[l][n]>999.0) continue;
            float dev=abs(eTablePair[k][l][m][n]-ak[m]-bl[n]);
            if(dev>maxdev){
              maxdev=dev;
            }
          }
        }

        if(maxdev<=threshold){
          adjMatrix[i][j]=adjMatrix[j][i]=0;
          for(int m=0;m<nrots[k];m++){
            if(eTableSelf[k][m]>999.) continue;
            eTableSelf[k][m]+=ak[m];
          }
          for(int n=0;n<nrots[l];n++){
            if(eTableSelf[l][n]>999.) continue;
            eTableSelf[l][n]+=bl[n];
          }
          //cout<<"edge between site "<<k<<" and "<<l<<" is decomposed with threshold "<<threshold<<endl;
          nedgeremoved++;
        }
      }
    }
  }
  cout<<"edge decomposition with threshold "<<threshold<<"... #edges removed: "<<nedgeremoved<<endl;

  int resremoved=0;
  for(int i=0;i<adjMatrix.size();i++){
    int k=unfixres[i];
    int nonzero=0;
    for(int j=0;j<adjMatrix.size();j++){
      if(adjMatrix[i][j]!=0){
        nonzero++;
        break;
      }
    }
    if(nonzero==0){
      int best;
      float emin=1e8;
      for(int j=0;j<nrots[k];j++){
        if(eTableSelf[k][j]>999.0) continue;
        if(eTableSelf[k][j]<emin){
          emin=eTableSelf[k][j];
          best=j;
        }
      }
      bestrot[k]=best;
      nrots[k]=1;
      //cout<<"residue "<<k<<" eliminated by edge decomposition"<<endl;
      resremoved++;
    }
  }
  cout<<"remove residues with no edge ... "<<"#residues removed: "<<resremoved<<endl;
}


void Solution::Search()
{
  IV1 pos;
  for(int i=0;i<nres;i++){
    if(nrots[i]>1) pos.push_back(i);
  }
  if(DEESearch(pos)==true){
    bool hardmode=false;
    float threshold=0.5;
LOOP:
    int nunfix=unfixres.size();
    IV1 visited(nunfix,0);
    IV2 adjMatrix(nunfix,visited);
    IV2 flagMatrix(nunfix,visited);
    //construct subgraphs
    ConstructAdjMatrix(nunfix,adjMatrix);
    if(hardmode==true){
      hardmode=false;
      GraphEdgeDecomposition(adjMatrix,threshold);
      threshold*=2;
    }
    ConstructSubgraphs(nunfix,visited,adjMatrix,flagMatrix);
    ShowGraphs();

    for(int i1=0; i1<graphs.size(); i1++){
      Graph &graph=graphs[i1];
      tree.Subgraph2TreeDecomposition(i1,graph);
      int treewidth=tree.CheckTreewidth();
      cout<<"#treewidth = "<<treewidth<<endl;
      if(treewidth<=TREEWIDTH_CUT){
        TreeDecompositionBottomToTopCalcEnergy();
        RootBagFindGMEC(tree.connBags[0]);
        for(set<int>::iterator it=tree.connBags[0].childBagIdx.begin();it!=tree.connBags[0].childBagIdx.end();it++){
          TreeDecompositionTopToBottomAssignRotamer(tree.connBags[0],tree.connBags[*it]);
        }
        for(Graph::iterator it2=graph.begin();it2!=graph.end();it2++){
          int site=it2->first;
          for(vector<int>::iterator it=unfixres.begin();it!=unfixres.end();it++){
            if(*it==site){
              it=unfixres.erase(it);
              break;
            }
          }
        }
        cout<<"current tree has been solved"<<endl;
      }
      else{
        hardmode=true;
        cout<<"current tree is hard, to be solved later"<<endl;
      }
      TreeDecompositionRelease();
      graph.clear();
    }
    graphs.clear();
    if(hardmode) goto LOOP;
  }

  //cout<<"optimal rotamer index:"<<endl;
  for(int i=0;i<nres;i++){
    if(nrots[i]==1){
      Pick(i,bestrot[i]);
      //cout<<i<<"\t"<<optRotIdx[i]<<endl;
    }
  }

  pdb=stru;
  stru.clear();
  sc.clear();
}
