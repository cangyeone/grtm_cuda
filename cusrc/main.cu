#include <stdio.h> 
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "popm/popm.h"
int main(int argc, char **argv){

    green("../cusrc/example/crust1.0.txt", "out4.txt", "signal.txt");
    return 0; 
}