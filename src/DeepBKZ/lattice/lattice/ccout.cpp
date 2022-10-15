#ifndef _inc_ccout_cpp
#define _inc_ccout_cpp

//for organized display
class dstreambuf : public streambuf {
    public:

    char buf[4096];
    
    dstreambuf() {
        char* pstart = buf;
        char* pend = buf + 4096;
        setp(pstart,pend);
        setg(pstart,pstart,pend);
    }
    
    virtual int sync() {
        *pptr() = '\0';
        int ilen = pptr() - pbase();

        pbump (pbase() - pptr());
        struct winsize w;
        ioctl(0,TIOCGWINSZ,&w);

        int cpos = 0;
        
        for (int i=0;i<ilen;i++) {
            if (buf[i] == '\r') {
                cpos = 0;
                cout << '\r';
                cout.flush();
            } else
            if (buf[i] == '\n') {
                cout << endl;
                cpos= 0;
            } else {
                if (cpos < w.ws_col) {
                    printf("%c",buf[i]);
                }
                cpos++;
            }
        }

        return 0;
    }
    
};
dstreambuf myds;
std::iostream ccout(&myds);

#endif