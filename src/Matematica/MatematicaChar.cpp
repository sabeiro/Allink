#include <Matematica.h>

uint sTable(char *String){
  if(strncmp(String,"s",1)){
    return 0x03C3;
  }
  else if(strncmp(String,"S",1)){
    return 0x03A3;
  }



}
uint sTable(char Char){
  return (uint) Char + 0x0350;
}
