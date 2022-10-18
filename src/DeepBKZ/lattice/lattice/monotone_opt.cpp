#ifndef _inc_monotone_opt
#define _inc_monotone_opt

//To find the inverse value of a monotone function
template <typename T> struct ValueTable {
    
        std::vector<std::vector< T> > data; 
        T direction = 1;    //+1=monotone inc. -1=monotone dec
        
        void add(T a,T b) {
            std::vector<T> temp;
            temp.resize(2);
            temp[0] = a;
            temp[1] = b*direction;
            data.push_back(temp);
            std::sort(data.begin(),data.end());
        }
        
        void cleartable() {
            data.clear();
        }
        
        void display() {
            for (int i=0;i<data.size();i++) {
                ccout << data[i][0] << " " << data[i][1] << endl;
            }
        }
        
        T nextcandidate(T target,T& range) {
            target *= direction;
            //assumed to be sorted
            int size = data.size();
            if (size==1) return data[0][0] * 2;
            range = data[size-1][0] - data[0][0];
            if (data[0][1] > target ) return data[0][0] * 0.5;
            if (data[size-1][1] < target) return data[size-1][0] * 2.0;
            
            int rsize = max(size /4,1);
            int i = size /2;
            while(1) {
                if (i<=0) i=0;
                if (i>=size-2) i = size-2;
                //ccout << "nextcandidate"  << i << " " << data[i][1] << "(" << data[i][0] << ") " << target << " " << data[i+1][1] << "(" << data[i+1][0] << ") " << endl;
                
                if ((data[i][1] <= target) && (target <= data[i+1][1] )) {
                    T ratio = (target-data[i][1]) / (data[i+1][1]-data[i][1]);
                    range = data[i+1][0]-data[i][0];
                    if (ratio < 0.1) ratio = 0.1; 
                    if (ratio > 0.9) ratio = 0.9; 
                    return data[i][0] + ratio * (data[i+1][0]-data[i][0]);
                }
                if (target <= data[i][1]) {
                    i-=rsize;
                    rsize = max(rsize/2,1);
                }                
                if (target > data[i+1][1]) {
                    i+=rsize;
                    rsize = max(rsize/2,1);
                }                
            }
        }
    
};

template <typename Tkey,typename Tval> struct COptpair {
    Tkey k;
    Tval v;
};

template <typename Tkey,typename Tval> bool operator <(const COptpair<Tkey,Tval>& a,const COptpair<Tkey,Tval>& b) {
    if (a.v < b.v) return true;
    if (a.v > b.v) return false;
    if (a.k < b.k) return true;
    return false;
}

template <typename Tkey,typename Tval> bool operator >(const COptpair<Tkey,Tval>& a,const COptpair<Tkey,Tval>& b) {
    if (a.v > b.v) return true;
    if (a.v < b.v) return false;
    if (a.k > b.k) return true;
    return false;
}

//Optimizing 
template <typename Tkey,typename Tval> struct COptTable {

    std::vector<COptpair<Tkey,Tval > > data; 
    int direction=+1;   //+1=maximize -1=minimize
    void add(Tkey k,Tval v) {
        data.push_back({k,v});
    }    
    
    bool valexist(Tkey k) {
        for (int i=0;i<data.size();i++) {
            if (data[i].k == k) return true;
        }
        return false;
    }
    
    Tkey minkey() {
        Tval minval = data[0].v;
        Tkey ret = data[0].k;
        for (int i=1;i<data.size();i++) {
            if (minval > data[i].v) {
                minval = data[i].v;
                ret = data[i].k;
            }
        }
        return ret;
    }

    Tval minvalue() {
        Tval minval = data[0].v;
        for (int i=1;i<data.size();i++) {
            if (minval > data[i].v) {
                minval = data[i].v;
            }
        }
        return minval;
    }

    Tkey nextvaluerandom(Tkey& range) {
        if (data.size()<=2) {
            ccout << "COptTable: error? data.size=" << data.size() << endl;
        }
        if (direction==-1) {
            std::sort(data.begin(),data.end());
        } else {
            std::sort(data.begin(),data.end(),std::greater<COptpair<Tkey,Tval> >());
        }

        //quadratic equation 
        Tkey x2 = data[0].k;
        Tkey x1 = -1;
        Tkey x3 = -1;
        for (int i=0;i<data.size();i++) {
            if (x2 < data[i].k) {
                if (x3<0) x3 = data[i].k;
                if (x3 >  data[i].k) x3 = data[i].k;
            }
            if (x2 > data[i].k) {
                if (x1<0) x1 = data[i].k;
                if (x1 <  data[i].k) x1 = data[i].k;
            }
        }
        Tkey xmax = max(x1,max(x2,x3));
        Tkey xmin = min(x1,min(x2,x3));
        Tkey xrnd;
        while (1) {
            xrnd = xmin + (1.0 * (rand()%10000)/10000.0)*(xmax-xmin);
            if (valexist(xrnd)==false) return xrnd;
        }
    }

    Tkey nextvalue(Tkey& range) {
        if (data.size()<=2) {
            ccout << "COptTable: error? data.size=" << data.size() << endl;
        }
        if (direction==-1) {
            std::sort(data.begin(),data.end());
        } else {
            std::sort(data.begin(),data.end(),std::greater<COptpair<Tkey,Tval> >());
        }

        //quadratic equation 
        Tkey x2 = data[0].k;
        Tval y2 = data[0].v;
        Tkey x1 = x2;
        Tkey x3 = x2;
        Tval y1,y3;
        for (int i=1;i<data.size();i++) {
            if (x2 < data[i].k) {
                if ((x3==x2) || (x3 >  data[i].k)) {
                    x3 = data[i].k;
                    y3 = data[i].v;
                }
            }
            if (x2 > data[i].k) {
                if ((x1==x2) || (x1  < data[i].k)) {
                    x1 = data[i].k;
                    y1 = data[i].v;
                }
            }
        }
        range = max(x1,max(x2,x3)) - min(x1,min(x2,x3));
        Tkey ret;
        if ( x2-x1 < x3-x2 ) {
            ret = (x2+x3)/2;
        } else {
            ret = (x1+x2)/2;
        }
        if (ret < min(x1,min(x2,x3))) ret = min(x1,min(x2,x3));
        if (ret > max(x1,max(x2,x3)) ) ret = max(x1,max(x2,x3)) ;
        return ret;
    }
    
    void display() {
        for (int i=0;i<data.size();i++) {
            ccout << data[i].k << " " << data[i].v << endl;
        }
    }
};


template <typename T1,typename T2> struct keypair {
    T1 key;
    T2 value;
};

template <typename T1,typename T2> bool operator <( const keypair<T1,T2>& a,const  keypair<T1,T2>& b ) {
    return (a.key<b.key);
    
}

template <typename T1,typename T2> class OptConvex {
    
    std::vector< keypair<T1,T2> > table;
    T1 range_l;
    T1 range_h;
    int internal_counter;
    
    public:
        
    OptConvex(T1 low,T1 high) {
        range_l = low;
        range_h = high;
        internal_counter = 0;
    }
    
    void add(T1 key,T2 value) {
        for (int i=0;i<table.size();i++) {
            if (table[i].key==key) return;
        }
        keypair<T1,T2> temp;
        temp.key = key;
        temp.value=value;
        table.push_back(temp);
    }
    
    T2 getmin() {
        T2 ret = table[0].value;
        for (int i=1;i<table.size();i++) {
            ret = min(ret,table[i].value);
        }
        return ret;
    }

    T1 getminkey() {
        T2 minvalue = table[0].value;
        T1 ret = table[0].key;
        for (int i=1;i<table.size();i++) {
            if (minvalue > table[i].value) {
                minvalue = table[i].value;
                ret = table[i].key;
            }
        }
        return ret;
    }

    
    T1 next(T1& range) {
        range = range_h-range_l;
        if (table.size()==0) return range_l;
        if (table.size()==1) return range_h;
        if (table.size()==2) return (range_h+range_l)/2;
        
        std::sort(table.begin(),table.end());
        int minindex=0;
        T2 minvalue = table[0].value;
        
        for (int i=0;i<table.size();i++) {
            if (minvalue > table[i].value) {
                minindex = i;
                minvalue = table[i].value;
            }
        }
 
        char direction = 1;
        if ((internal_counter++)%2==0) {   
            direction = -1;
        }
        if (minindex==0) {
            direction = 1;
        }
        if (minindex==table.size()-1) {
            direction = -1;
        }
        range = abs(table[minindex].key - table[minindex + direction].key); 
        return (table[minindex].key + table[minindex + direction].key)/2;
    }
    
    
};




#endif