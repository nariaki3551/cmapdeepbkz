#ifndef _inc_time_tools
#define _inc_time_tools

//get wallclock time in the format "yyyy/mm/dd/hh:mm:ss"
std::string timestr() {
    std::ostringstream ret;
    time_t timer;
    timer = time(NULL);
    struct tm *lct;
    lct = localtime(&timer);
    ret << lct->tm_year + 1900 << "/";
    if (lct->tm_mon+1<10) ret << "0";
    ret << lct->tm_mon+1 << "/";
    if (lct->tm_mday<10) ret << "0";
    ret << lct->tm_mday << "/";
    if (lct->tm_hour<10) ret << "0";
    ret << lct->tm_hour << ":";
    if (lct->tm_min<10) ret << "0";
    ret << lct->tm_min << ":";
    if (lct->tm_sec<10) ret << "0";
    ret << lct->tm_sec;
    return ret.str();
}

//timer 
double gettimeofday_sec() { 
    struct timeval tv; 
    gettimeofday(&tv, NULL); 
    return tv.tv_sec + 0.000001 * tv.tv_usec; 
}

double getcputime() {
    return (double)clock()/(double)CLOCKS_PER_SEC;
}



#endif