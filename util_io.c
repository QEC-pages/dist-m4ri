#include "util_io.h"
// #include "util.h"

//extern params_t pio;
params_io_t pio={1,3,0,0};

int local_io_init(int argc, char **argv, csr_t **spaP, csr_t **spaG){
  int i;
  char *finP, *finG;
  //  memset(&prm,0,sizeof(params_t));
  int dbg=0;
  int swit=0;
  finP="tryX.mtx";
  finG="tryZ.mtx";

  for (i=1; i<argc;i++) /* scan arguments for help message */
    if((strcmp(argv[i],"--help")==0)||(strcmp(argv[i],"-h")==0))
      return i; 

  for(i=1; i<argc; i++){
    if(sscanf(argv[i],"debug=%d",& dbg)==1){
      if(dbg==0)
	pio.debug=0;
      else{
	pio.debug|=dbg;
	printf("# read %s, debug=%d octal=%o\n",argv[i],pio.debug,pio.debug);
      }
    }					
    else if (sscanf(argv[i],"css=%d",&dbg)==1){
      pio.css=dbg;
      if (pio.debug)
	printf("# read %s, css=%d\n",argv[i],pio.css);
    }
    else if (sscanf(argv[i],"switch=%d",&dbg)==1){ /* switch P and G files */
      swit=dbg;
      if (pio.debug)
	printf("# read %s, switch=%d\n",argv[i],swit);
    }
    else if (0==strncmp(argv[i],"finP=",5)){
      if (swit){
	finG=argv[i]+5; 
	if (pio.debug)
	  printf("# read %s; since switch=%d: finG=%s\n",argv[i],swit,finG);
      }
      else{ 
	finP=argv[i]+5; 
	if (pio.debug)
	  printf("# read %s; finP=%s\n",argv[i],finP);
      }
    }
    else if (0==strncmp(argv[i],"finG=",5)){
      if (swit){
	finP=argv[i]+5; 
	if (pio.debug)
	  printf("# read %s; since switch=%d: finP=%s\n",argv[i],swit,finP);
      }
      else{ 
	finG=argv[i]+5; 
	if (pio.debug)
	  printf("# read %s; finG=%s\n",argv[i],finG);
      }
    }
    else if (0==strncmp(argv[i],"fin=",4)){
      int len = strlen(argv[i]+4);
      char *s = (char *) malloc((len+6)*sizeof(char));
      sprintf(s,"%s%s",argv[i]+4,swit>0?"Z.mtx":"X.mtx");
      finG=s;
      s = (char *) malloc((len+6)*sizeof(char));
      sprintf(s,"%s%s",argv[i]+4,swit>0?"X.mtx":"Z.mtx");
      finP=s;
      if (pio.debug)
	printf("# read %s; since switch=%d assigning \n# finG=%s\n# finP=%s\n",
	       argv[i],swit,finG,finP);
    }
    else{ /* unrecognized option */
      if (pio.debug)
	printf("# end of prm_IO at position %d: \"%s\" \n",i,argv[i]);
      break; /* this is the end of the arguments that can be interpreted here */ 
    }
      
  }

  //  printf("####### here i=%d ###\n",i);
  //  if((strcmp(argv[i],"--help")!=0)&&(strcmp(argv[i],"-h")!=0)){
    /* read matrices unless help message is requested */

  if (*spaP!=NULL) *spaP=csr_free(*spaP);
  *spaP=csr_mm_read(finP,NULL,0); 
  if (*spaG!=NULL) *spaG=csr_free(*spaP);
  *spaG=csr_mm_read(finG,NULL,0);
  rci_t n = (*spaP)-> cols;
  if (n != (*spaG) -> cols)
    ERROR("Column count mismatch in P and G matrices: %d != %d",
	  (*spaP)-> cols, (*spaG)->cols);
  pio.n=n; 
  pio.n0=n;
  if (pio.css!=1)
    ERROR("Non-CSS codes are currently not supported, css=%d",pio.css);
  
  return i;
}

