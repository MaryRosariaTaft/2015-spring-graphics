#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

int main(){
  
  int fd;
  char buf[32];
  int xres, yres, max_val;

  xres = 1920/2;
  yres = 1080/2;
  max_val = 256;
  
  fd = open("pic.ppm", O_CREAT | O_TRUNC | O_WRONLY, 0664);
  sprintf(buf, "P3\n%d %d\n%d\n", xres, yres, max_val);
  write(fd, buf, strlen(buf));
  
  int i, j; //counters
  int r, g, b; //colors
  
  for(i=0; i<yres; i++){
    for(j=0; j<xres; j++){

      if(i<(yres/10) || i>(9*yres/10)){
	r = 123;
	g = 0;
	b = j%256;
      }else{
	r = (i*i)%max_val;
	g = (j*j)%max_val;
	b = (i+j)%max_val;
      }

      sprintf(buf, "%d %d %d\t", r, g, b);
      write(fd, buf, strlen(buf));

    }
    
  }

  close(fd);

  return 0;
  
}
