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
#include "Search.h"
#include <sstream>

using namespace std;
string PROGRAM_PATH=(string)".";
string ROTLIB2010=(string)"dun2010bbdep.bin";

int main(int argc,char** argv)
{
  cout<<"###########################################################################"<<endl;
  cout<<"                    FASPR (Version 20200309)                 "<<endl;
  cout<<"  A method for fast and accurate protein side-chain packing, "<<endl;
  cout<<"which is an important problem in protein structure prediction"<<endl;
  cout<<"and protein design."<<endl;
  cout<<endl;
  cout<<"Copyright (c) 2020 Xiaoqiang Huang"<<endl;
  cout<<"Yang Zhang Lab"<<endl;
  cout<<"Dept. of Computational Medicine and Bioinformatics"<<endl;
  cout<<"Medical School"<<endl;
  cout<<"University of Michigan"<<endl;
  cout<<"Email:tommyhuangthu@foxmail.com, xiaoqiah@umich.edu"<<endl;
  cout<<"###########################################################################"<<endl;
  clock_t start,finish;
  float duration;
  start = clock();

  if(argc<2){
    cout<<"Usage: ./FASPR -i input.pdb -o output.pdb\n";
    cout<<"[-s sequence.txt] to load a sequence file\n";
    return 0;
  }

  /**********************************************************/
  /* get program path and check if rotamer library exists   */
  /**********************************************************/
  char fullpath[2048];
  strcpy(fullpath,argv[0]);
  for(int i = strlen(fullpath); i >= 0; --i){
    if(fullpath[i]=='/' || fullpath[i]=='\\'){
      fullpath[i + 1]='\0';
      break;
    }
  }
  PROGRAM_PATH=(string)fullpath;
  string rotfile=PROGRAM_PATH+"/"+ROTLIB2010;
  fstream infile(rotfile.c_str(),ios::in|ios::binary);
  if(!infile){
    cerr<<"error! cannot find rotamer library "<<ROTLIB2010<<endl;
    exit(0);
  }
  else{
    infile.close();
  }

  string pdbin=(string)"example/1mol.pdb";
  string pdbout=(string)"example/1mol_FASPR.pdb";
  string seqfile=(string)"void";

  bool sflag=false;
  int i;
  for(i=1;i<argc-1;i++){
   if(argv[i][0]=='-'){
     if(argv[i][1]=='i'){
       i++;
       pdbin=argv[i];
     }
     else if(argv[i][1]=='o'){
       i++;
       pdbout=argv[i];
     }
     else if(argv[i][1]=='s'){
       i++;
       seqfile=argv[i];
       sflag=true;
     }
   }
  }
  
  Solution faspr;
  faspr.ReadPDB(pdbin);
  if(sflag) faspr.LoadSeq(seqfile);
  else faspr.LoadSeq();
  faspr.BuildSidechain();
  faspr.CalcSelfEnergy();
  faspr.CalcPairEnergy();
  faspr.Search();
  faspr.WritePDB(pdbout);
  finish = clock();
  duration = (float)(finish-start)/CLOCKS_PER_SEC;
  cout<<"#computational time: "<<duration<<" seconds"<<endl;

  return 0;
}
