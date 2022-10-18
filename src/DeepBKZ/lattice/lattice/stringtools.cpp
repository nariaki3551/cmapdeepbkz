#ifndef _inc_string_tools
#define _inc_string_tools

template <typename T> std::string to_stdstring(T x) {
    stringstream ss;
    ss << x;
    return ss.str();
}
  
template <typename T> std::string to_hexstring(T x) {
    stringstream ss;
    ss << std::hex << x;
    return ss.str();
}


void ExtractOptions(std::map<std::string,std::string>& options,std::string opts) {
    
    std::vector<std::string> args;
    std::vector<std::string> args2;
    boost::algorithm::split(args,opts,boost::is_any_of(" "));
    
    for (int i=0;i<args.size();i++) {
        //ccout << args[i] << endl;
        boost::algorithm::split(args2,args[i],boost::is_any_of("="));
        std::string aa;
        if (args2.size()>=2) {
            aa = args2[1];
            for (int j=2;j<args2.size();j++) {
                aa += "=" + args2[j];
            }
        } else {
            aa = "true";
        }
        options[args2[0]] = aa;
    }
}


#endif