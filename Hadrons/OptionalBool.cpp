struct InvalidOptionalException : public std::exception
{
	const char * what () const throw ()
    {
    	return "Access of uninitialised value";
    }
};

class optionalBool {
private:
    bool value_;
    bool valid_ = false;
public:
    optionalBool() : valid_(false) {}
    optionalBool(bool value) : value_(value) , valid_(true) { 
    }

    bool operator= (bool value) {
        value_ = value;
        valid_ = true;
        return this->value_;
    }

    bool value() const {
        if (valid_)
            return value_;

        throw InvalidOptionalException();
    }

    operator bool() const {
        return value();
    }

    friend std::ostream& operator<<(std::ostream& os, const optionalBool& b) {
        os << b.value();
        return os;
    }

    friend std::istream& operator>>(std::istream& os, optionalBool& b) {
        bool value;
        os >> value;
        b = value;
        return os;
    }
};