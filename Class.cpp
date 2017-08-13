#include "MagneticField.h"

int main()
{
  MagneticField * obj = new MagneticField(10., 10., 10., 5., 5., 5., 0.1, 10., 1000);
  obj->CreateRootOutput("field_output");
  
  return 0;
}
