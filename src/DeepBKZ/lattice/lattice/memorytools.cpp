#ifndef _inc_memory_tools_cpp
#define _inc_memory_tools_cpp

//Subroutines for shared memory for progressive bkz library
namespace shared_memory {

    struct typeandid {
        unsigned int thash;      //hash(typeid.name())
        int id;
    };
    bool operator < (const typeandid& a,const typeandid& b) {
        if (a.thash < b.thash) return true;
        if (a.thash == b.thash) {
             if (a.id < b.id) return true;
        }
        return false;
    }

    struct darray {
        void* d=NULL;
        int size1=0;
        int size2=0;
        int size3=0;
    };

    std::map<typeandid,darray> memory;
    
    template <typename T> void* allocate1(int id,int size,int vl=0) {
        //set vl!=0 when debug
        struct typeandid ti;
        std::hash<std::string> shash;
        std::string tname = typeid(T).name();
#ifdef __no__multithreads
        ti.thash = shash(tname);
#else
        ti.thash = shash(tname) + omp_get_thread_num();
#endif        
        ti.id = id;
        if (vl!=0) {
            ccout << "---------------------------" << endl;
            ccout << "memoryalloc: type=" << typeid(T).name() << " id=" << id << endl;
            ccout << "thash=" << ti.thash << endl;
        }        
        darray dd;
        if (memory.find(ti) == memory.end()) {
            if (vl!=0) ccout << "not found" << endl;
            //not exist
        } else {
            dd = memory.at(ti);
        }
        if (vl!=0) ccout << "allocated size=" << dd.size1 << endl;
        
        if (dd.size1 < size) {
            if (dd.d!=NULL) delete [] (T*)dd.d;
            dd.d = new T[size];
            T* pointer = (T*)dd.d;
            for (int i=0;i<size;i++) (*pointer++) = 0;
            dd.size1 = size;
            if (vl!=0) ccout << "allocate new memory of type=" << typeid(T).name() << " size=" << size << endl;
        }
        memory[ti] = dd;
        return (void*)dd.d;
    }

    template <typename T> void* allocate1_init(int id,int size,int vl,T init_val) {
        //set vl!=0 when debug
        struct typeandid ti;
        std::hash<std::string> shash;
        std::string tname = typeid(T).name();
#ifdef __no__multithreads
        ti.thash = shash(tname);
#else
        ti.thash = shash(tname) + omp_get_thread_num();
#endif        
        ti.id = id;
        if (vl!=0) {
            ccout << "---------------------------" << endl;
            ccout << "memoryalloc: type=" << typeid(T).name() << " id=" << id << endl;
            ccout << "thash=" << ti.thash << endl;
        }        
        darray dd;
        if (memory.find(ti) == memory.end()) {
            if (vl!=0) ccout << "not found" << endl;
            //not exist
        } else {
            dd = memory.at(ti);
        }
        if (vl!=0) ccout << "allocated size=" << dd.size1 << endl;
        
        if (dd.size1 < size) {
            if (dd.d!=NULL) delete [] (T*)dd.d;
            dd.d = new T[size];
            T* pointer = (T*)dd.d;
            for (int i=0;i<size;i++) (*pointer++) = 0;
            dd.size1 = size;
            if (vl!=0) ccout << "allocate new memory of type=" << typeid(T).name() << " size=" << size << endl;
        }
        memory[ti] = dd;
        return (void*)dd.d;
    }

    template <typename T> void* allocate2(int id,int s1,int s2,int vl=0) {
        struct typeandid ti;
        std::hash<std::string> shash;
        std::string tname = typeid(T).name();
#ifdef __no__multithreads
        ti.thash = shash(tname + "dim2");
#else
        ti.thash = shash(tname + "dim2") + omp_get_thread_num();;
#endif        
        ti.id = id;
        if (vl!=0) {
            ccout << "---------------------------" << endl;
            ccout << "memoryalloc: type=" << typeid(T).name() << " id=" << id << endl;
            ccout << "thash=" << ti.thash << endl;
        }        
        darray dd;
        if (memory.find(ti) == memory.end()) {
            if (vl!=0) ccout << "not found" << endl;
            //not exist
        } else {
            dd = memory.at(ti);
        }
        if (vl!=0) ccout << "allocated size=" << dd.size1  << endl;
        
        if ((dd.size1 < s1) || (dd.size2 < s2)){
            if (dd.d!=NULL) {
                T** temp = (T**)dd.d;
                for (int i=0;i<dd.size1;i++) delete [] temp[i];
                delete [] temp;
            }
            T** temp;
            temp = new T*[s1];
            for (int i=0;i<s1;i++) {
                temp[i] = new T[s2];
            }
            
            dd.d = (void*)temp;
            dd.size1 = s1;
            dd.size2 = s2;
            if (vl!=0) ccout << "allocate new memory of type=" << typeid(T).name() << " size=" << s1 << ","  << s2<< endl;
        }
        memory[ti] = dd;
        return (void*)dd.d;
    }

    template <typename T> void* allocate2_init(int id,int s1,int s2,int vl=0) {
        struct typeandid ti;
        std::hash<std::string> shash;
        std::string tname = typeid(T).name();
#ifdef __no__multithreads
        ti.thash = shash(tname + "dim2");
#else
        ti.thash = shash(tname + "dim2") + omp_get_thread_num();;
#endif        
        ti.id = id;
        if (vl!=0) {
            ccout << "---------------------------" << endl;
            ccout << "memoryalloc: type=" << typeid(T).name() << " id=" << id << endl;
            ccout << "thash=" << ti.thash << endl;
        }        
        darray dd;
        if (memory.find(ti) == memory.end()) {
            if (vl!=0) ccout << "not found" << endl;
            //not exist
        } else {
            dd = memory.at(ti);
        }
        if (vl!=0) ccout << "allocated size=" << dd.size1  << endl;
        
        if ((dd.size1 < s1) || (dd.size2 < s2)){
            if (dd.d!=NULL) {
                T** temp = (T**)dd.d;
                for (int i=0;i<dd.size1;i++) delete [] temp[i];
                delete [] temp;
            }
            T** temp;
            temp = new T*[s1];
            for (int i=0;i<s1;i++) {
                temp[i] = new T[s2];
                for (int j=0;j<s2;j++) {
                    temp[i][j] = 0;
                }
            }
            
            dd.d = (void*)temp;
            dd.size1 = s1;
            dd.size2 = s2;
            if (vl!=0) ccout << "allocate new memory of type=" << typeid(T).name() << " size=" << s1 << ","  << s2<< endl;
        }
        memory[ti] = dd;
        return (void*)dd.d;
    }


    template <typename T> void* allocate3(int id,int s1,int s2,int s3,int vl=0) {
        struct typeandid ti;
        std::hash<std::string> shash;
        std::string tname = typeid(T).name();
#ifdef __no__multithreads
        ti.thash = shash(tname + "dim3");
#else
        ti.thash = shash(tname + "dim3") + omp_get_thread_num();    
#endif        
        ti.id = id;
        if (vl!=0) {
            ccout << "---------------------------" << endl;
            ccout << "memoryalloc: type=" << typeid(T).name() << " id=" << id << endl;
            ccout << "thash=" << ti.thash << endl;
        }        
        darray dd;
        if (memory.find(ti) == memory.end()) {
            if (vl!=0) ccout << "not found" << endl;
            //not exist
        } else {
            dd = memory.at(ti);
        }
        if (vl!=0) ccout << "allocated size=" << dd.size1 << endl;
        
        if ((dd.size1 < s1) || (dd.size2 < s2) || (dd.size3 < s3)) {
            if (dd.d!=NULL) {
                T*** temp = (T***)dd.d;
                for (int i=0;i<dd.size1;i++) {
                    for (int j=0;j<dd.size2;j++) {
                        delete [] temp[i][j];
                    }
                    delete [] temp[i];
                }
                delete [] temp;
            }
            T*** temp;
            s1 -= (s1%8);     s1 += 8;
            s2 -= (s1%8);     s2 += 8;
            s3 -= (s1%8);     s3 += 8;

            temp = new T**[s1];
            for (int i=0;i<s1;i++) {
                temp[i] = new T*[s2];
                for (int j=0;j<s2;j++) {
                    temp[i][j] = new T[s3];
                }
            }
            dd.d = (void*)temp;
            dd.size1 = s1;
            dd.size2 = s2;
            dd.size3 = s3;
            if (vl!=0) ccout << "allocate new memory of type=" << typeid(T).name() << " size=" << s1 << "," << s2 << "," << s3 << endl;
        }
        memory[ti] = dd;
        return (void*)dd.d;
    }

    template <typename T> void  release(int id,int size,int vl=0) {
        //Todo
    }
    
}



#endif