
/**
 * \file gff.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "regexp.h"
#include "gff.h"

using namespace fsa;

const std::string GFF::comment_char = "#";
const std::string GFF::undef_char = ".";

const std::string GFF::attributes_split_char = ";";
const std::string GFF::attributes_assign_char = "=";
const std::string GFF::attributes_list_char = ",";

const std::string GFF::key_id = "ID";
const std::string GFF::key_name = "Name";

Regexp GFF::re_key_value = "^([^=]+)=(.*)$";

bool ciCharLess(char c1, char c2)
{
    return tolower(static_cast<unsigned char>(c1)) <
           tolower(static_cast<unsigned char>(c2));
}

bool ciStringCompare::operator() (const std::string& s1, const std::string& s2) const
{
    return std::lexicographical_compare(s1.begin(),
                                    s1.end(),
                                    s2.begin(),
                                    s2.end(),
                                    ciCharLess);
}

void GFF::from_string (const std::string& str) {

    std::vector<std::string> tokens = Util::tokenize_string (str, "\t"); // hold tab-separated tokens

    // now parse tokens

    // sanity checks
    if (tokens.size() != 9) {
        cerr << "Not a GFF-format line: " << str << endl;
		return;
    }

    seqid = tokens[0];
    source = tokens[1];
    type = tokens[2];
    start = static_cast<unsigned> (atoi (tokens[3].c_str()));
    end = static_cast<unsigned> (atoi (tokens[4].c_str()));
    score = (tokens[5] != "" && tokens[5] != GFF::undef_char) ? (float)atof (tokens[5].c_str()) : -1.f;
    strand = (tokens[6])[0];
    phase = (tokens[7] != "" && tokens[7] != GFF::undef_char) ? static_cast<unsigned> (atoi (tokens[7].c_str())) : 3;

    parse_attributes_string (tokens[8]);

}

void GFF::parse_attributes_string (const std::string& str) {

    std::vector<std::string> key_value_pairs = Util::tokenize_string (str, GFF::attributes_split_char);
    for (std::vector<std::string>::const_iterator key_value = key_value_pairs.begin(); key_value != key_value_pairs.end(); ++key_value) {
        if (re_key_value.Match (key_value->c_str())) {
            std::string key = Util::trim(re_key_value[1]);
            std::vector<std::string> values = Util::tokenize_string (re_key_value[2], GFF::attributes_list_char);
            if (values.size())
                attributes_ordering.push_back (key);
            for (std::vector<std::string>::const_iterator value = values.begin(); value != values.end(); ++value)
                attributes[key].push_back (Util::trim(*value));
        }
    }

}

std::string GFF::get_attributes_string() const {

    std::string s;

    // iterate over keys
    for (std::vector<std::string>::const_iterator key = attributes_ordering.begin(); key != attributes_ordering.end(); ++key) {
        std::map<std::string, std::vector<std::string> >::const_iterator key_value = attributes.find (*key);
        // iterate over values
        if (key_value->second.size())
            s += *key + GFF::attributes_assign_char;
        for (std::vector<std::string>::const_iterator value = key_value->second.begin(); value != key_value->second.end(); ++value)
            s += *value + GFF::attributes_list_char;
        if (s.length())
            s.erase (s.end() - 1, s.end()); // strip off final comma
        s += GFF::attributes_split_char;
    }

    if (!s.length())
        s += GFF::undef_char;

    return s;

}

std::string GFF::to_string() const {

    std::stringstream ss;

    ss << seqid << "\t"
       << (source != "" ? source : GFF::undef_char) << "\t"
       << (type != "" ? type : GFF::undef_char) << "\t"
       << start << "\t" << end << "\t";
    if (score != -1.) { ss << score; }
    else { ss << GFF::undef_char; }
    ss << "\t";
    ss << strand << "\t";
    if (phase != 3) { ss << phase; }
    else { ss << GFF::undef_char; }
    ss << "\t";
    ss << get_attributes_string();

    return ss.str();

}

void GFF_database::from_file (const std::string& filename) {

    Regexp re_gff ("^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*$");

    // open file
    std::ifstream filestream;
    filestream.open (filename.c_str(), std::ios::in);
    if (!filestream.is_open()) {
        cerr << "Couldn't open file '" << filename << "' for reading." << endl;
        exit (1);
    }

    cerr << "Reading GFF file '" << filename << "'...";

    // parse file
    std::string line;
    while (!filestream.eof()) {
        getline (filestream, line);
		
		std::string::size_type cr_ret = line.rfind('\r');
		if (cr_ret != std::string::npos)
			line.resize(cr_ret);
		
        Util::chomp (line);
        // skip lines which don't seem to be in GFF format
        if (!line.length() || !re_gff.Match (line.c_str()))
            continue;
        GFF gff;
		
        gff.from_string (line);
		if (gff.seqid == "")
			continue;
        store_entry (gff);
        __maxlen = (gff.end - gff.start + 1 > __maxlen) ? gff.end - gff.start + 1 : __maxlen;
    }

    cerr << "done." << endl;

    // clean up
    filestream.close();

}

void GFF_database::write (std::ostream& o) const {

    for (std::vector<GFF>::const_iterator gff = this->begin(); gff != this->end(); ++gff)
        o << *gff;

}

void GFF_database::append (const GFF_database& gff_db) {

    for (std::vector<GFF>::const_iterator gff = gff_db.begin(); gff != gff_db.end(); ++gff)
        this->store_entry (*gff);
    __maxlen = (gff_db.maxlen() > __maxlen) ? gff_db.maxlen() : __maxlen;

}

void GFF_database::create_unique_ids() {

    if (!size())
        return;

    size_t width = static_cast<size_t> (log10 (size()));
    std::string format = "%" + Util::to_string (width) + "d";
    for (size_t i = 0; i < size(); ++i) {
        std::string id (width - static_cast<size_t> (log10 (i)), '0'); // hack to format the IDs nicely
        id += Util::to_string (i);
        __entries[i].set_id (id);
    }

}

void GFF_database::sort_entries() {

    std::sort (__entries.begin(), __entries.end(), GFF::GFF_less());
    __is_sorted = true;

}

GFF_database GFF_database::chromosome_features (const std::string& chromosome) const {

    if (!__is_sorted) {
        cerr << "You must call sort_entries before attepting use this function." << endl;
        exit (1);
    }

    GFF_database features;

    // check for empty interval or no features
    if (!size())
        return features;

    // initialize dummy object for comparison
    GFF g;
    g.seqid = chromosome;
    //Util::strip_leading_chr (g.seqid);

    g.start = 1; // convert to 1-based coordinates
    g.end = 1;
    std::vector<GFF>::const_iterator gff_lower = lower_bound (__entries.begin(),
                                                              __entries.end(),
                                                              g,
                                                              GFF::GFF_less());

    g.start = 0xFFFFFFFE + 1; // convert to 1-based coordinates
    g.end =   0xFFFFFFFE + 1;
    std::vector<GFF>::const_iterator gff_upper = upper_bound (__entries.begin(),
                                                              __entries.end(),
                                                              g,
                                                              GFF::GFF_less());

    // now iterate backwards through the sorted entries,
    // checking the second condition (1.end >= 2.start)
    // as we go
    while (gff_lower < gff_upper) {

        if (gff_lower->seqid != chromosome)
            continue;

        features.store_entry (*gff_lower);

        ++gff_lower;
    }

    // sort nicely
    features.sort_entries();

    return features;

}


GFF_database GFF_database::intersect_genomic_interval (const std::string& chromosome,
                                                       const unsigned& start, const unsigned& end) const {

    if (!__is_sorted) {
        cerr << "You must call sort_entries before attepting use this function." << endl;
        exit (1);
    }

    GFF_database intersections;

    // check for empty interval or no features
    if (end < start || !size())
        return intersections;

    // initialize dummy object for comparison
    GFF g;
    g.seqid = chromosome;
    Util::strip_leading_chr (g.seqid);

    // two intervals 1 and 2 intersect iff
    // 1.start <= 2.end && 1.end >= 2.start
    // features in database are 1; query interval is 2

    // check first condition (1.start <= 2.end)
    g.start = end + 1; // convert to 1-based coordinates
    g.end = end + 1;
    std::vector<GFF>::const_iterator gff_upper = upper_bound (__entries.begin(),
                                                              __entries.end(),
                                                              g,
                                                              GFF::GFF_less());

    // decrement element to try to get to the upper-bound entry
    // which (possibly) intersects the query interval
    if (gff_upper != __entries.begin())
        --gff_upper;

    // now check that we're on the correct chromosome
    if (gff_upper->seqid != chromosome)
        return intersections;

    // now iterate backwards through the sorted entries,
    // checking the second condition (1.end >= 2.start)
    // as we go
    while (gff_upper != __entries.begin() - 1) {

        // if we're no longer on the correct chromosome, then we're done
        if (gff_upper->seqid != chromosome)
            break;

        if (gff_upper->start > end) {
            cerr << "query is " << chromosome << " " << start << ", " << end << endl;
            cerr << "gff_upper - 2 = " << *(gff_upper - 2);
            cerr << "gff_upper - 1 = " << *(gff_upper - 1);
            cerr << "gff_upper = " << *gff_upper;
            cerr << "gff_upper + 1 = " << *(gff_upper + 1);
            cerr << "gff_upper + 2 = " << *(gff_upper + 2);
        }
        assert (gff_upper->start <= end);

        // does it intersect the query interval?
        if (gff_upper->end >= start)
            intersections.store_entry (*gff_upper);

        // check whether we're done
        // (i.e., whether given __maxlen, there cannot be any more features close
        // enough to the query interval to intersect it)
        if (gff_upper->start + __maxlen < start)
            break;

        --gff_upper;
    }

    // sort nicely
    intersections.sort_entries();

    return intersections;

}
