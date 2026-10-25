

/*
Contains the STRING_REF struct
*/

#pragma once

#ifndef COMMON_HPP
#define COMMON_HPP 

#include <string>
#include <cstring>
#include <cstdint>
#include <iostream>



using namespace std;


struct STRING_REF{
    STRING_REF(const string*s,int st,int len):str(s),start(st),length(len){}
    STRING_REF() = default;
    const string* str;
    uint16_t start;
    uint16_t length;
    bool operator==(const STRING_REF&other)const{
        if(length != other.length){return false;}
        if(str == other.str && start == other.start){return true;}
        return(
            std::memcmp(str->data() + start, other.str->data() + other.start, length) == 0
        );
    }
    const char& operator[](size_t i)const{
        return(str->data()[start+i]);
    }
    void print_last_char()const{
        cout.write(str->data()+start+length-1,1);
    }
    friend ostream& operator<<(ostream&out,const STRING_REF& ref);
};


ostream& operator<<(ostream&out,const STRING_REF& ref){
    if(ref.str && ref.length>0){
        out.write(ref.str->data() + ref.start, ref.length);
    }
    return out;
}



#endif