
/**
 * \file gff.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef SEQ_GFF_INCLUDED
#define SEQ_GFF_INCLUDED

#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <stdint.h>

#include "misc.h"

namespace fsa {

/* This is from Meyers, "Effective C++" 3rd ed. */
struct ciStringCompare
{
    bool operator() (const std::string& s1, const std::string& s2) const;
};

/**
 * \brief Bare-bones representation of a single GFF entry.
 *
 * This class repeatedly assumes that the seqid field holds the
 * chromosome which the feature is on.
 * Note that GFF coordinates are always 1-based and fully-closed.
 * This class stores coordinates for GFF features accordingly
 * (in contrast to the 0-based indexing used throughout other
 * parts of this code).
 * See http://www.sequenceontology.org/gff3.shtml.
 */
struct GFF {

public:

    typedef std::map<std::string, std::vector<std::string>, ciStringCompare > AttributeTable;
    std::string seqid;        ///< field 0
    std::string source;       ///< field 1
    std::string type;
    unsigned start;
    unsigned end;
    float score;
    char strand;
    unsigned phase;
    std::vector<std::string> attributes_ordering;                     ///< sorted list of keys in attributes field
    AttributeTable attributes;  ///< map of key, value pairs

    /**
     * \brief Default constructor.
     */
    GFF()
        : seqid (""), source (""), type (""),
        score (-1.f), strand ('.'), phase (3)
    {
    }

    /**
     * \brief Constructor.
     *
     * Initialize from a scored interval.
     */
    GFF (std::string seqid, unsigned start, unsigned end, float score)
        : seqid (seqid), source (""), type (""),
        start (start), end (end),
        score (score), strand ('.'), phase (3)
    {
    }

    /**
     * \brief Add a (possibly-new) key, value pair to the attributes field.
     *
     * If a value for the key already exists, then the new value is appended.
     */
    template<typename K, typename V>
    inline void add_value (const K& key, const V& value);

    /**
     * \brief Set a (possibly-new) key, value pair to the attributes field.
     *
     * If a value for the key already exists, then the old value is cleared
     * and replaced with the new value.
     */
    template<typename K, typename V>
    inline void set_value (const K& key, const V& value);

    /**
     * \brief Set name in attributes field.
     */
    inline void set_name (const std::string& name);

    /**
     * \brief Set ID in attributes field.
     */
    inline void set_id (const std::string& id);

    /**
     * \brief Initialize from a GFF-formatted std::string.
     */
    void from_string (const std::string& str);

    /**
     * \brief Write to GFF-formatted string.
     */
    std::string to_string() const;

    /**
     * \brief Output operator.
     *
     * Prints undef_char for undefined fields.
     */
    friend std::ostream& operator<< (std::ostream& o, const GFF &gff) {
        o << gff.to_string() << endl;
        return o;
    }


    static const std::string comment_char;                    ///< comment character
    static const std::string undef_char;                      ///< character for undefined fields

    static const std::string attributes_split_char;           ///< split key, value pairs in attributes field
    static const std::string attributes_assign_char;          ///< assign value to key in attributes field
    static const std::string attributes_list_char;            ///< separate values for a single key in attributes field

    static const std::string key_id;                          ///< GFF key for ID
    static const std::string key_name;                        ///< GFF key for Name

    /**
     * \brief Function object for binary comparison of GFF objects (sort by seqid, start, end).
     */
    struct GFF_less : std::binary_function<GFF, GFF, bool> {
public:
        bool operator() (const GFF &l, const GFF &r) const {
            if (l.seqid == r.seqid) {
                if (l.start == r.start)
                    return l.end < r.end;
                return l.start < r.start;
            }
            return l.seqid < r.seqid;
        }
    };


protected:

    /**
     * \brief Get string for the attributes field.
     */
    std::string get_attributes_string() const;

private:

    /**
     * \brief Set start coordinate.
     */
    inline void set_start (const unsigned s);

    /**
     * \brief Parse attributes string into key, value pairs.
     */
    void parse_attributes_string (const std::string& str);

    static Regexp re_key_value;                                  ///< match key, value pairs in attributes field

};

/**
 * \brief Represent a GFF file.
 */
struct GFF_database {

    typedef std::vector<GFF>::iterator iterator;
    typedef std::vector<GFF>::const_iterator const_iterator;
    /**
     * \brief Constructor.
     */
    GFF_database()
        : __maxlen (0) {
    }

    /**
     * \brief Load from a file.
     */
    void from_file (const std::string& filename);

    /**
     * \brief Write.
     */
    void write (std::ostream& o) const;

    /**
     * \brief Store an entry.
     */
    inline void store_entry (const GFF& gff);

    /**
     * \brief Append a GFF_database to this one.
     */
    void append (const GFF_database& gff_db);

    /**
     * \brief Create nicely-formatted numerical IDs for all entries.
     */
    void create_unique_ids();

    /**
     * \brief Sort entries in database.
     *
     * This must be called in order for intersect_genomic_interval to work properly.
     */
    void sort_entries();

    /**
     * \brief Find GFF features which intersect the passed interval.
     *
     * Assumes that the seqid field holds the chromosome
     * and that the passed interval coordinates are 0-based
     * and fully closed.  Note the contrast with the 1-based
     * coordinates used for the GFF features themselves.
     * Note that the entries MUST be sorted in order for this to function properly.
     * \see GFF_database::sort_entries
     */
    GFF_database intersect_genomic_interval (const std::string& chromosome,
                                             const unsigned& start, const unsigned& end) const;

    GFF_database chromosome_features (const std::string& chromosome) const;

    /**
     * \brief Maximum length of entry.
     */
    size_t maxlen() const {
        return __maxlen;
    }

    /**
     * \brief Number of entries.
     */
    size_t size() const {
        return __entries.size();
    }

    /**
     * \brief Data access operator.
     */
    const GFF& operator[] (const size_t &i) const {
        assert (i < __entries.size());
        return __entries[i];
    }


    /**
     * \brief Get iterator to start of __entries.
     */
    iterator begin() {
        return __entries.begin();
    }

    /**
     * \brief Get const_iterator to start of __entries.
     */
    const_iterator begin() const {
        return __entries.begin();
    }

    /**
     * \brief Get iterator to end of __entries.
     */
    iterator end() {
        return __entries.end();
    }

    /**
     * \brief Get const_iterator to end of __entries.
     */
    const_iterator end() const {
        return __entries.end();
    }

private:

    std::vector<GFF> __entries; ///< individual GFF entries in database

    size_t __maxlen;            ///< maximum length of an entry in the database
    bool __is_sorted;           ///< have the entries been sorted?

};



/****************************************
* Function definitions.
****************************************/



inline void GFF::set_start (const unsigned s) {
    if (s >= 1)
        start = s;
    else {
        cerr << "Setting start coordinate " << s << " to 0." << endl;
        start = 0;
    }
}

template<typename K, typename V>
inline void GFF::add_value (const K& key, const V& value) {

    // use stringstream to convert key and value to string
    std::string key_str = Util::to_string (key);
    std::string value_str = Util::to_string (value);

    // now store
    // if we already have this key, then record new key
    if (attributes.find (key_str) == attributes.end())
        attributes_ordering.push_back (key_str);
    // append value
    attributes[key_str].push_back (value_str);

}

template<typename K, typename V>
inline void GFF::set_value (const K& key, const V& value) {

    // use stringstream to convert key and value to string
    std::string key_str = Util::to_string (key);
    std::string value_str = Util::to_string (value);

    // now store
    // if we already have this key, then clear old value and record new value
    AttributeTable::iterator key_itr = attributes.find(key_str);

    if (key_itr != attributes.end()) {
        key_itr->second.clear();
        key_itr->second.push_back (value_str);
    }
    else {
        attributes_ordering.push_back (key_str);

        // TODO: insert hint here?
        attributes[key_str].push_back (value_str);
    }

}

inline void GFF::set_name (const std::string& name) {
    set_value (key_name, name);
}

inline void GFF::set_id (const std::string& id) {
    set_value (key_id, id);
}

inline void GFF_database::store_entry (const GFF& gff) {

    __entries.push_back (gff);
    __maxlen = (gff.end - gff.start + 1 > __maxlen) ? gff.end - gff.start + 1 : __maxlen;

}

}

#endif /* SEQ_GFF_INCLUDED */
