#ifndef UTILS_H
#define UTILS_H

// static c-style: these functions' scope are limited to file in which 
// #include is used

static int fcount_pts(const char *infile) {

  //const char *c_str = infile.c_str();
  //FILE *fvec = fopen(c_str, "r");
  std::ifstream fvec(infile);
  // counts the number of lines in the given file not beginning with "#"
  // to get the number of points
  if(fvec==NULL) { 
    fprintf(stderr,"Could not open file %s",infile);
    exit(1);
  }
  std::string line;
  int points = 0;
  while (std::getline(fvec, line))
  {  
      if (line.compare(0, 1, "#") != 0)
          ++points;
  }
  //fclose(fvec);

  return points;
}

static void FindMax2(const double *errors,
                 const int len,
                 double &worst_err,
                 int &worst_app)
{
  worst_err = 0.0;
  for(int i = 0; i < len; i++) {
    if(worst_err < errors[i]) {
      worst_err = errors[i];
      worst_app = i;
    }
  }
}





#endif

