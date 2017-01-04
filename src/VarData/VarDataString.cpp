/***********************************************************************
VarData: This  Program reads and writes a specific file format 
storing all the information relative to a set of equal structure
polymers to the CHAIN, PART and GENERAL structures. It provides 
two different ways to backfold the coordinates, a fanction that 
creates an initial system with different option and some function
for the data analisys. The first calculate the distribution of the
monomer in the box, the second the distribution of the bonds.
Copyright (C) 2008 by Giovanni Marelli <sabeiro@virgilio.it>


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
***********************************************************************/
#include "../include/VarData.h"
#include <stdarg.h>

bool VarData::ReadString(const char *String,char *cLine,double *Value){
  if(String == NULL) return 1;
  int iLen = (int) (strlen(cLine));
  char cVar[STRSIZE];
  char cVal[STRSIZE];
  int Paren[2];
  for(int i=0;i<iLen;i++){
    if(cLine[i]=='#')//comment
      break;
    if(cLine[i]==String[0]){
      Paren[0] = i;
      for(int j=i;j<iLen;j++){
	if(cLine[j]==' '){
	  Paren[1] = j;
	  memset(cVar,0,STRSIZE);
	  strncpy(cVar,cLine+Paren[0],Paren[1]-Paren[0]);
	  //printf("%s %s\n",cVar,String);
	  for(int q =0;q<Paren[1]-Paren[0];q++) cVar[q] = cLine[q+Paren[0]];
	  if( (strcmp(cVar,String))==0 ){
	    for(int l=j;l<iLen;l++){
	      if(cLine[l]!=' '){
		Paren[0] = l;
		for(int h=l;h<iLen;h++){
		  if(cLine[h]==' ' || cLine[h]=='\n'){
		    Paren[1] = h;
		    //sprintf(cVal,"");
		    strncpy(cVal,cLine+Paren[0],Paren[1]-Paren[0]);
		    //printf("\t%d %d %s\n",Paren[0],Paren[1],cVal); 
		    sprintf(cVal+Paren[1]-Paren[0],"   ");
		    //strncpy(cVal+Paren[1]-Paren[0],cSpa,3);
		    sscanf(cVal,"%lf",Value);
		    //printf("%s Vale %f\n",cVar,*Value);
		    return 1;
		  }
		}
		break;
	      }
	    }
	    break;
	  }
	  else break;
	}
      }
      //break;
    }
  }
  return 0;
}
bool VarData::ReadString(const char *String,double *Value,char *Line){
  //printf("--%s\n",Line);
  char *pLine = Line;
  char cVar[STRSIZE];
  char cVal[STRSIZE];
  int Paren[2];
  if(pLine!=NULL){
    int iLen = (int) (strlen(Line));
    //printf("%s\n",Line);
    printf("%s->%s",String,Line);
    for(int i=0;i<iLen;i++){
      if(pLine[i]=='#')//comment
	break;
       if(pLine[i]==String[0]){
	Paren[0] = i;
	for(int j=i;j<iLen;j++){
	  if(pLine[j]==' '){
	    Paren[1] = j;
	    memset(cVar,STRSIZE,sizeof(char));
	    strncpy(cVar,Line+Paren[0],Paren[1]-Paren[0]);
	    //for(int q =0;q<Paren[1]-Paren[0];q++) cVar[q] = Line[q+Paren[0]];
	    printf("  %d %d %s   %d\n",Paren[0],Paren[1],cVar,strlen(cVar)); 
	    if( (strcmp(cVar,String))==0 ){
	      //	      printf("%s\n",Line);
	      //printf("  %d %d %s\n",Paren[0],Paren[1],cVar); 
	      for(int l=j;l<iLen;l++){
		if(pLine[l]!=' '){
		  Paren[0] = l;
		  for(int h=l;h<iLen;h++){
		    if(pLine[h]==' ' || pLine[h]=='\n'){
		      Paren[1] = h;
		      //sprintf(cVal,"");
		      strncpy(cVal,Line+Paren[0],Paren[1]-Paren[0]);
		      //printf("\t%d %d %s\n",Paren[0],Paren[1],cVal); 
		      sprintf(cVal+Paren[1]-Paren[0],"   ");
		      sscanf(cVal,"%lf",Value);
		      printf("%s Vale %f\n",cVar,*Value);
		      return 1;
		    }
		  }
		  break;
		}
	      }
	      break;
	    }
	    else break;
	  }
	}
	//break;
      }
    }
  }
  return 0;
}
int VarData::ReadVal(char *pLine,double *Value){
  if(pLine==NULL) return 0;
  char cVar[STRSIZE];
  char cVal[20];
  int Paren[2];
  int Incr = 0;
  *Value = 0.;
  int iLen = (int) (strlen(pLine));
  memset(cVar,0,STRSIZE);
  for(int i=0;i<iLen;i++){
    if(pLine[i]=='#')//comment
      return 0;
    if(pLine[i]!=' '){
      Paren[0] = i;
      for(int j=i;j<iLen;j++){
	if(pLine[j]==' ' || pLine[j]=='\n'){
	  Incr = Paren[1] = j;
	  strncpy(cVar,pLine+Paren[0],Paren[1]-Paren[0]);
	  sprintf(cVal+Paren[1]-Paren[0],"   ");
	  sscanf(cVar,"%lf",Value);
	  //printf("  %d %s %lf\n",Incr,cVar,*Value);
	  pLine += Incr;
	  return Incr;
	}
      }
    }
  }
  return 0;
}
int VarData::Fetch(char *str,char *mask,int NArg,double *Val){
  int sLen = strlen(str);
  int sPos = 0;
  int IfContinue = 0;
  // if(!strncmp(str,mask,strlen(mask))) return 0;
  // char *pOpen = strpbrk(str,"([");
  // char *pClose = strpbrk(str,")]");
  // if(pOpen == NULL || pClose == NULL) return 0;
  // sPos = pOpen - str +1;
  // sLen = pClose - str;
  for(int s=0;s<sLen;s++){
    if(!strncmp(str+s,mask,strlen(mask)) ){
      //if(str[s] == mask[0]){
      sPos = s;
      for(int ss=s;ss<sLen;ss++){
  	if(str[ss] == '(' || str[ss] == '['){
  	  sPos = ss+1;
  	  for(int sss=ss;sss<sLen;sss++){
  	    if(str[sss] == ')' || str[sss] == ']'){
  	      sLen = sss;
  	      IfContinue = 1;
  	      break;
  	    }
  	  }
  	  break;
  	}
      }
      break;
    }
  }
  if(!IfContinue) return 0;
  //char *sNumber = (char*)calloc(160,sizeof(char));
  //strncpy(sNumber,str+sPos,sLen-sPos);
  char *sNumber = str+sPos;
  int sNum = 0;
  for(int a=0;a<NArg;a++){
    sscanf(sNumber+sNum,"%lf",Val+a);
    //FIXME: two blank spaces/one letter
    for(int s=sNum+1;s<sLen-sPos;s++){
      if(sNumber[s] == ' '){
	for(int ss=s;ss<sLen-sPos;ss++){
	  if(sNumber[ss] != ' '){
	    sNum = ss;
	    break;
	  }
	}
	break;
      }
    }
    //printf("%s-%d) %s vale %lf in %d\n",mask,a,sNumber+sNum,*(Val+a),sPos);
  }
  //free(sNumber);
  return sPos;
}
int VarData::BraketPos(char *str,char *mask,int *sPos,int *sLen){
  if(!strncmp(str,mask,strlen(mask))) return 1;
  char *pInit = strpbrk(str,mask);
  if(pInit == NULL) return 1;
  int InitPos = pInit - str + 1;
  char *pOpen = strpbrk(str+InitPos,"([{");
  char *pClose = strpbrk(str+InitPos,"}])");
  //printf("%s %d %d %d %s\n",mask,pInit,*sPos,*sLen,str+InitPos);
  if(pOpen == NULL || pClose == NULL) return 1;
  *sPos = pOpen - str +1;
  *sLen = pClose - str;
  //printf("%s %d %d %d %s\n",mask,pInit,*sPos,*sLen,str+InitPos);
  return 0;
}
int VarData::Fetch(char *str,char *mask,char *fmt, ... ){
  int sLen = 0;
  int sPos = 0;
  int IfContinue = BraketPos(str,mask,&sPos,&sLen);
  if(IfContinue) return 0;
  char Field[120];
  memset(Field,' ',120*sizeof(char));
  //char *Field = (char *)calloc(sLen-sPos+1,sizeof(char));
  strncpy(Field,str+sPos,sLen-sPos);
  va_list args;
  va_start(args,fmt);
  vsscanf(Field,fmt,args);
  //  printf("Sotto %s trovato %s come %s\n",mask,Field,fmt);
  va_end(args);
  //free(Field);
  return sPos;
}
// //FIXME: fa seg fault
// int VarData::Fetch(char *str,char *mask,char *fmt, ... ){
//   int sLen = 0;
//   int sPos = 0;
//   int IfContinue = BraketPos(str,mask,&sPos,&sLen);
//   if(IfContinue) return 0;
//   char *Field = (char *)calloc(sLen-sPos+1,sizeof(char));
//   if(Field == NULL){
//     printf("Fetch non allocated in %s %s %s\n",mask,str,fmt);
//     return 0;
//   }
//   // char Field[60];
//   // memset(Field,' ',60*sizeof(char));
//   strncpy(Field,str+sPos,sLen-sPos);
//   int ret = 0;
//   va_list args;
//   va_start(args,fmt);
//   ret = vsscanf(Field,fmt,args);
//   //printf("Sotto %s trovato %s come %s\n",mask,Field,fmt);
//   va_end(args);
//   free(Field);
//   return sPos;
// }
