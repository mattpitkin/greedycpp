#ifndef UTILS_H
#define UTILS_H

static int fcount_pts(const char *infile) {
    //const char *c_str = infile.c_str();
    //FILE *fvec = fopen(c_str, "r");
    std::ifstream fvec(infile);
    // counts the number of lines in the given file not beginning with "#"
    // to get the number of points
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


#endif
