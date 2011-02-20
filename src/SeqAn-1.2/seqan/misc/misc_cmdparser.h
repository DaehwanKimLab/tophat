 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_MISC_CMDPARSER
#define SEQAN_MISC_CMDPARSER

#include <sstream>
#include <seqan/map.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//  TODO:
//      * support some more formating options
//      * store/return error code (invalid argument, invalid option, etc.)
//      * support named arguments (e.g. <ARG1> -> <INPUT FILE>)
//////////////////////////////////////////////////////////////////////////////
template<typename TChar>
inline bool
_isDigit(TChar const c)
{
    return (c >= '0') && (c <= '9');
}

template<typename TString>
inline bool
_isDouble(TString const s)
{
    bool _dot = true;
    unsigned l = length(s);
    unsigned i = 0;

    // skip leading sign
    if(s[i] == '-') ++i;
    while(i < l){
        if(!_isDigit(s[i])){
            if(s[i] == '.' && _dot){
                _dot = false;
            }else return false;
        }
        ++i;
    }
    return true;
}

template<typename TString>
inline bool
_isInt(TString const s)
{
    unsigned l = length(s);
    unsigned i = 0;
    // skip leading sign
    if (s[i] == '-') ++i;
    while(i < l){
        if(!_isDigit(s[i])) return false;
        ++i;
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////

struct OptionType{
    enum {
        Boolean = 1,		// option needs no argument, value is true iff given on command line
        String = 2,			// argument is a string
        Int = 4,			// ... an integer
        Double = 8,			// ... a float
        Mandatory = 16,		// option must be set
		Label = 32,			// automatically print a label for the argument(s) on the help screen
		List = 64,			// option is a list of values
        Hidden = 128		// hide this option from the help screen
	};
};

//////////////////////////////////////////////////////////////////////////////

/**
.Class.CommandLineOption:
..cat:Miscellaneous
..summary:Stores information for a specific Commandline Option
..signature:CommandLineOption
*/
class CommandLineOption
{
public:
    CharString			longName;			// long option name
    CharString			shortName;			// short option name
	CharString			arguments;			// argument names seperated by spaces

    CharString			helpText;			// option description
    int					optionType;			// option type
	int					argumentsPerOption;	// number of arguments per option
	
	String<CharString>	defaultValue;
	String<CharString>	value;

    CommandLineOption() {}

    CommandLineOption(
		CharString const & _short,
		CharString const & _long,
		CharString const & _help,
		int _type
	) :
		longName(_long),
		shortName(_short),
		helpText(_help),
		optionType(_type),
		argumentsPerOption(1)
	{
	}

    CommandLineOption(
		CharString const & _short,
		CharString const & _long,
		int _argumentsPerOption,
		CharString const & _help,
		int _type
	) :
		longName(_long),
		shortName(_short),
		helpText(_help),
		optionType(_type),
		argumentsPerOption(_argumentsPerOption)
	{
	}

	template <typename TValue>
    CommandLineOption(
		CharString const & _short,
		CharString const & _long,
		int _argumentsPerOption,
		CharString const & _help,
		int _type,
		TValue const & _default
	) :
		longName(_long),
		shortName(_short),
		helpText(_help),
		optionType(_type),
		argumentsPerOption(_argumentsPerOption)
	{
		std::stringstream strm;
		strm << _default;
		appendValue(defaultValue, strm.str());
		append(helpText, " (default ");
		append(helpText, strm.str());
		appendValue(helpText, ')');
	}

	template <typename TValue>
    CommandLineOption(
		CharString const & _short,
		CharString const & _long,
		CharString const & _help,
		int _type,
		TValue const & _default
	) :
		longName(_long),
		shortName(_short),
		helpText(_help),
		optionType(_type),
		argumentsPerOption(1)
	{
		std::stringstream strm;
		strm << _default;
		appendValue(defaultValue, strm.str());
		append(helpText, " (default ");
		append(helpText, strm.str());
		appendValue(helpText, ')');
	}

/**.Memfunc.CommandLineOption#CommandLineOption:
..class:Class.CommandLineOption
..summary:Constructor
..signature:CommandLineOption ()
..signature:CommandLineOption (shortName,longName,helpText,type)
..param.shortName:A @Shortcut.CharString@ containing the short option identifier (e.g. $"h"$ for the $-h/--help$ option).
...remarks:Note that the leading "-" is not passed.
..param.longName:A @Shortcut.CharString@ containing the long option identifier (e.g. $"help"$ for the $-h/--help$ option).
...type:Shortcut.CharString
...remarks:Note that the leading "--" is not passed.
..param.helpText:A @Shortcut.CharString@ containing the help text associated with this option.
...type:Shortcut.CharString
*/
};

//////////////////////////////////////////////////////////////////////////////

inline CommandLineOption
addArgumentText(CommandLineOption const & opt, CharString const & text)
{
	CommandLineOption temp = opt;
	temp.arguments = " ";
	append(temp.arguments, text);
	return temp;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.longName:
..summary:Returns the long option name of a @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:longName(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:A @Shortcut.CharString@ holding the long name of the CommandLine Option (e.g. $help$ in case of $-h/--help$)
..remarks:The result type is @Shortcut.CharString@.
*/
inline CharString &
longName(CommandLineOption & me){
    return me.longName;
}

inline const CharString &
longName(CommandLineOption const & me){
    return me.longName;
}

/**
.Function.setLongName:
..summary:Sets the long option name of a @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:setLongName(option,newName)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.newName:A @Shortcut.CharString@ containing the new long name of the option.
...type:Shortcut.CharString
*/
inline void
setLongName(CommandLineOption & me, CharString const & newName){
    me.longName = newName;
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.shortName:
..summary:Returns the short option name of a @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:shortName(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:A @Shortcut.CharString@ holding the short name of the CommandLine Option (e.g. $h$ in case of $-h/--help$)
..remarks:The result type is @Shortcut.CharString@.
*/
inline CharString &
shortName(CommandLineOption & me){
    return me.shortName;
}

inline const CharString &
shortName(CommandLineOption const & me){
    return me.shortName;
}

/**
.Function.setShortName:
..summary:Sets the short option name of a @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:setShortName(option,newName)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.newName:A @Shortcut.CharString@ containing the new short name of the option.
*/
inline void
setShortName(CommandLineOption & me, CharString const & newName){
    me.shortName = newName;
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.helpText:
..summary:Returns the help text associated with the @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:helpText(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:A @Shortcut.CharString@ holding the help text of the CommandLine Option
..remarks:The result type is @Shortcut.CharString@.
*/
inline CharString &
helpText(CommandLineOption & me){
    return me.helpText;
}

inline const CharString &
helpText(CommandLineOption const & me){
    return me.helpText;
}

/**
.Function.setHelpText:
..summary:Sets the help text associated with the @Class.CommandLineOption@ object
..cat:Miscellaneous
..signature:setHelpText(option,newHelpText)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.newHelpText:A @Shortcut.CharString@ containing the new help text.
...type:Shortcut.CharString
*/
inline void
setHelpText(CommandLineOption & me, CharString const & newHelp){
    me.helpText = newHelp;
}

//////////////////////////////////////////////////////////////////////////////

inline bool
isStringOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::String) != 0;
}

inline bool
isBooleanOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Boolean) != 0;
}

inline bool
isDoubleOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Double) != 0;
}

inline bool
isIntOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Int) != 0;
}

inline bool
isHiddenOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Hidden) != 0;
}

inline bool
isOptionMandatory(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Mandatory) != 0;
}

inline bool
isLabelOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Label) != 0;
}

inline bool
isOptionList(CommandLineOption const & me)
{
    return (me.optionType & OptionType::List) != 0;
}

inline void
setOptionType(CommandLineOption & me, const int _newOpt)
{
    me.optionType = _newOpt;
}

//////////////////////////////////////////////////////////////////////////////

inline CharString
argumentText(CommandLineOption const & me)
{	
	if (empty(me.arguments))
	{
		CharString label;
		if (isLabelOption(me))
		{
			if (isStringOption(me))
				label = " STR";
			else if (isIntOption(me) || isDoubleOption(me))
				label = " NUM";
			
			if (me.argumentsPerOption >= 2)
			{
				std::stringstream strm;
				if (!empty(label))
					for (int i = 0; i < me.argumentsPerOption; ++i)
						strm << label << (i + 1);
				return strm.str();
			}
		}
		return label;
    }
	else
		return me.arguments;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream>
inline void
_writeOptName(TStream & target, CommandLineOption const & me)
{
    _streamWrite(target, empty(shortName(me)) ? "" : "-");
    _streamWrite(target, shortName(me));
    _streamWrite(target, (empty(shortName(me)) || empty(longName(me))) ? "" : ", ");
    if (!empty(longName(me)))
    {
        _streamWrite(target, "--");
        _streamWrite(target, longName(me));
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream>
inline void
write(TStream & target, CommandLineOption const & me)
{
    _streamPut(target,'\t');
    _writeOptName(target, me);
    _streamPut(target,'\t');
    _streamPut(target,'\t');
    _streamWrite(target,me.helpText);
}

template <typename TStream>
inline TStream &
operator << (TStream & target, CommandLineOption const & source)
{
    write(target, source);
    return target;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Class.CommandLineParser:
..cat:Miscellaneous
..summary:Stores multiple @Class.CommandLineOption@ objects and parses the command line arguments for these options
..signature:CommandLineParser
*/
class CommandLineParser
{
public:
    typedef String<CommandLineOption>           TOptionMap;
    typedef Size<TOptionMap>::Type              TSize;
    
    typedef ::std::map<CharString, TSize>       TStringMap;
    typedef String<CharString>                  TValueMap;

    TStringMap           shortNameMap;
    TStringMap           longNameMap;
    TOptionMap           optionMap;
    
    unsigned             required_arguments;
    String<CharString>   arguments;
    CharString           appName;
	String<CharString>   titleText;
    String<CharString>   usageText;
	String<CharString>   versionText;

    unsigned line_width;
    unsigned padding_left;
	unsigned short_width;
	unsigned long_width;
	unsigned full_width;

	const CharString			null;			// empty return values
	const String<CharString>	nullSet;
	
/**.Memfunc.CommandLineParser#CommandLineParser:
..class:Class.CommandLineParser
..summary:Constructor
..signature:CommandLineParser ()
..signature:CommandLineParser (applicationName)
..param.applicationName:A @Shortcut.CharString@ containing the name of the application.
..remarks:If the name of the application is not passed to the constructor it will be extracted from the command line.
*/

    CommandLineParser()
        : required_arguments(0)
	{
		CommandLineOption opt("h","help","displays this help message",OptionType::Boolean);
		appendValue(optionMap,opt);
		insert(shortNameMap,"h",0);
		insert(longNameMap,"help",0);
		//insert(_long2ShortMap,longName(opt),shortName(opt));

		line_width   = 32;
		padding_left = 2;
		short_width  = 0;
		long_width   = 0;
		full_width   = 0;

		appName = "";
	}

    CommandLineParser(CharString appName)
        : required_arguments(0),appName(appName)
	{
		CommandLineOption opt("h","help","displays this help message",OptionType::Boolean);
		appendValue(optionMap,opt);
		insert(shortNameMap,"h",0);
		insert(longNameMap,"help",0);

		line_width   = 32;
		padding_left = 2;
		short_width  = 0;
		long_width   = 0;
		full_width   = 0;
	}
};

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addOption:
..summary:adds an instance of @Class.CommandLineOption@ to the @Class.CommandLineParser@
..cat:Miscellaneous
..signature:addOption(parser,option)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.option:The new @Class.CommandLineOption@ object that should be added.
...type:Class.CommandLineOption
*/
inline void
addOption(CommandLineParser & me, CommandLineOption const & opt)
{
	unsigned labelLen = length(argumentText(opt));
    appendValue(me.optionMap, opt);
    if (!empty(shortName(opt)))
	{
		insert(me.shortNameMap,shortName(opt), length(me.optionMap) - 1);
		unsigned width = 3 + length(shortName(opt));
		if (me.short_width < width)
			me.short_width = width;
		if (empty(longName(opt)))
		{
			width += 1 + length(argumentText(opt));
			if (me.full_width < width)
				me.full_width = width;
		}
	}
    if (!empty(longName(opt)))
	{
		insert(me.longNameMap,longName(opt), length(me.optionMap) - 1);
		unsigned width = 3 + length(longName(opt)) + labelLen;
		if (me.long_width < width)
			me.long_width = width;
	}
}

template <typename TString>
inline void
addLine(CommandLineParser & me, TString const & line)
{
	addOption(me, CommandLineOption("", "", line, 0));
}

template <typename TString>
inline void
addHelpLine(CommandLineParser & me, TString const & line)
{
	addOption(me, CommandLineOption("", "", line, 1));
}

template <typename TString>
inline void
addSection(CommandLineParser & me, TString const & line)
{
	addLine(me, "");
	addLine(me, line);
}

template <typename TString>
inline void
addTitleLine(CommandLineParser & me, TString const & line)
{
	appendValue(me.titleText, line);
}

template <typename TString>
inline void
addVersionLine(CommandLineParser & me, TString const & line)
{
	if (empty(me.versionText))
		addOption(me, CommandLineOption("V", "version", "print version information", OptionType::Boolean));
	appendValue(me.versionText, line);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.appendCmdLine:
..summary:adds a line of text to the help output of the @Class.CommandLineParser@
..cat:Miscellaneous
..signature:appendCmdLine(parser,text)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.text:A text line that will be added to the help output.
...type:Shortcut.CharString
*/
inline void
addUsageLine(CommandLineParser & me, CharString const & line)
{
    appendValue(me.usageText, line);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.hasOptionLong:
..summary:Returns true if the there is an option registered in the parser, that has the passed optionIdentifier
..cat:Miscellaneous
..signature:hasOptionLong(parser, optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the long option.
*/
inline bool 
hasOptionLong(CommandLineParser const & me, CharString const & _long)
{
    return hasKey(me.longNameMap, _long);
}


/**
.Function.hasOptionShort:
..summary:Returns true if the there is an option registered in the parser, that has the passed optionIdentifier
..cat:Miscellaneous
..signature:hasOptionShort(parser, optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the short option.
*/
inline bool 
hasOptionShort(CommandLineParser const & me, CharString const & _short)
{
    return hasKey(me.shortNameMap, _short);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.requiredArguments:
..summary:using this option you can define how many non parameterized options are required by your program.
..cat:Miscellaneous
..signature:requireRemainder(parser, count)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.count:A $unsigned int$ defining the amount of non-parameterized options requried by your program.
*/
inline void
requiredArguments(CommandLineParser & me, unsigned count)
{
    me.required_arguments = count;
}	


//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TStream>
inline void
_printStringSet(TStringSet const & set, TStream & target)
{
    for(unsigned r = 0; r < length(set); ++r)
    {
        _streamWrite(target, set[r]);
		_streamPut(target, '\n');
    }
}

template <typename TStream>
inline void
_usage(CommandLineParser const & me, TStream & target)
{
	_streamWrite(target, "Usage: ");
	if (empty(me.usageText))
	{
		_streamWrite(target, me.appName);
		_streamWrite(target, " [OPTION]... ");
		for (unsigned r = 0; r < me.required_arguments; ++r)
		{
			_streamWrite(target, "<ARG");
			_streamPutInt(target, r + 1);
			_streamWrite(target,"> ");
		}
		_streamPut(target,'\n');
	}
	else
	{
		for (unsigned r = 0; r < length(me.usageText); ++r)
		{
			if (r) _streamWrite(target, "       ");
			_streamWrite(target, me.appName);
			_streamPut(target, ' ');
			_streamWrite(target, me.usageText[r]);
			_streamPut(target,'\n');
		}
	}
}

template <typename TStream>
inline void
_title(CommandLineParser const & me, TStream & target)
{
	_printStringSet(me.titleText, target);
}

/**
.Function.shortHelp:
..summary:Prints a short help message for your parser to the stream
..cat:Miscellaneous
..signature:shortHelp(parser,stream)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
*/
template <typename TStream>
inline void
shortHelp(CommandLineParser const & me, TStream & target)
{
	_title(me, target);
    _usage(me, target);
    _streamWrite(target, "Try '");
    _streamWrite(target, me.appName);
    _streamWrite(target, " --help' for more information.\n");
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.help:
..summary:Prints the complete help message for your parser to the stream.
..cat:Miscellaneous
..signature:help(parser[,stream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
...default: $std::cerr$
*/
template <typename TStream>
inline void
help(CommandLineParser const & me, TStream & target)
{
	_title(me, target);
    _streamPut(target, '\n');
    _usage(me,target);
    _streamPut(target, '\n');

    for (unsigned o = 0; o < length(me.optionMap); ++o)
    {
        const CommandLineOption & opt = me.optionMap[o];
        if (isHiddenOption(opt)) continue;								// do not print hidden options
		
		if (opt.optionType > 0)
		{       
			unsigned s = 0;
			for (; s < me.padding_left; ++s)
				_streamPut(target, ' ');
			
			unsigned t1 = s + me.short_width;							// first tab
			unsigned t2 = _max(t1 + me.long_width, me.full_width) + 1;	// second tab (one extra space looks better)

			if (!empty(shortName(opt)))
			{
				_streamPut(target, '-');
				_streamWrite(target, shortName(opt));
				s += 1 + length(shortName(opt));
				if (!empty(longName(opt)))
				{
					_streamPut(target, ',');
					++s;
				} else {
					_streamWrite(target, argumentText(opt));
					s += length(argumentText(opt));
				}
			}
			
			for (; s < t1; ++s)
				_streamPut(target, ' ');
			
			if (!empty(longName(opt)))
			{
				_streamWrite(target, "--");
				_streamWrite(target, longName(opt));
				_streamWrite(target, argumentText(opt));
				s += 2 + length(longName(opt)) + length(argumentText(opt));
			}

			for (; s < t2; ++s)
				_streamPut(target, ' ');
		}

		_streamWrite(target, helpText(opt));

/*
        if (s < me.line_width){
			for (; s < me.line_width; ++s)
				_streamPut(target, ' ');
            _streamWrite(target, helpText(opt));
        }
        else
        {
            _streamPut(target, '\n');
            s = 0;
			for (; s < me.line_width; ++s)
				_streamPut(target, ' ');
            _streamWrite(target, helpText(opt));
        }
*/
        _streamPut(target, '\n');
    }
	_streamPut(target, '\n');
}

inline void
help(CommandLineParser const & me)
{
    help(me,::std::cerr);
}

/**
.Function.version:
..summary:Prints a version text to the stream.
..cat:Miscellaneous
..signature:version(parser[,stream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
...default: $std::cerr$
*/
template <typename TStream>
inline void
version(CommandLineParser const & me, TStream & target)
{
	_printStringSet(me.versionText, target);
}

inline void
version(CommandLineParser const & me)
{
    version(me,::std::cerr);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.isSetShort:
..summary:Returns true if the option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSetShort(parser,optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the short option.
*/
inline bool
isSetShort(CommandLineParser & me, CharString const & shortName)
{
    if (!hasKey(me.shortNameMap, shortName))
		return false; // this option does not exist

	// if value != "" -> value was set
	return !empty(me.optionMap[cargo(me.shortNameMap, shortName)].value);
}

/**
.Function.isSetLong:
..summary:Returns true if the option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSetLong(parser,optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the long option.
*/
inline bool
isSetLong(CommandLineParser & me,CharString const & longName)
{
    if (!hasKey(me.longNameMap, longName))
		return false; // this option does not exist

	// if value != "" -> value was set
	return !empty(me.optionMap[cargo(me.longNameMap, longName)].value);
}

//////////////////////////////////////////////////////////////////////////////

inline bool
_allMandatorySet(CommandLineParser const & me)
{
    for (unsigned o = 0; o < length(me.optionMap); ++o)
        if (empty(me.optionMap[o].value) && isOptionMandatory(me.optionMap[o])) return false;
    return true;
}

//////////////////////////////////////////////////////////////////////////////

inline CharString
_parseAppName(CharString const & candidate)
{
    int i = length(candidate) - 1;
	
    for(; i >= 0; --i)
        if (candidate[i] == '\\' || candidate[i] == '/') 
            break;

    return suffix(candidate, i + 1);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TErrorStream>
bool
_assignOptionValue(CommandLineParser & me, unsigned option_index, CharString const & val, unsigned argNo, TErrorStream & estream)
{
    // get the option object
    CommandLineOption & opt = me.optionMap[option_index];
    if(isDoubleOption(opt)){
        if(!_isDouble(val))
        {
            _streamWrite(estream,me.appName);
            _streamWrite(estream,": ");
            _streamWrite(estream, "\"");
            _streamWrite(estream, val);
            _streamWrite(estream, "\" is not a valid double value for '");
            _writeOptName(estream, opt);
            _streamWrite(estream, "'\n");
            return false;
        }
    }else if(isIntOption(opt)){
        if(!_isInt(val))
        {
            _streamWrite(estream,me.appName);
            _streamWrite(estream,": ");
            _streamWrite(estream, "\"");
            _streamWrite(estream, val);
            _streamWrite(estream, "\" is not a valid integer value for '");
            _writeOptName(estream, opt);
            _streamWrite(estream, "'\n");
            return false;
        }
    }
	if (isOptionList(opt))
		appendValue(opt.value, val, Generous());
	else
	{
		if (argNo == 0) clear(opt.value);
		appendValue(opt.value, val, Exact());
	}
    return true;
}

template<typename TErrorStream>
inline bool 
_assignOptionValue(CommandLineParser & me, unsigned option_index, CharString const & val, TErrorStream & estream)
{
	return _assignOptionValue(me, option_index, val, 0, estream);
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.parse:
..summary:Returns true if the option was set on the parsed command line.
..cat:Miscellaneous
..signature:parse(parser,argc,argv[,errorStream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.argc:Count of the objects on the command line.
..param.argv:Array of the different command line arguments ($const char *argv[]$). 
..param.errorStream:A stream where error messages are send too.
*/
template<typename TErrorStream>
bool
parse(CommandLineParser & me, int argc, const char *argv[], TErrorStream & estream)
{
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;
    // if the appName wasn't set .. parse from command line
    if (empty(me.appName)) me.appName = _parseAppName(argv[0]);

    for (int i = 1; i < argc; ++i) 
    {
        if (argv[i][0] == '-')  // this is possibly an option value
        {
            CharString inParam = argv[i];
            unsigned len = length(inParam);
            
            if (len == 1)
            {
                _streamWrite(estream,me.appName);
                _streamWrite(estream,": invalid option '-'\n");
                return false;
            }
            else if (inParam[1] != '-') // maybe a combination of multiple bool opts
            {
                for (unsigned s = 1; s < len; ++s)
				{
					unsigned e = len;
					for (; s < e; --e)
					{
						if (hasOptionShort(me, infix(inParam, s, e)))
						{
							TOptionPosition option_index = cargo(me.shortNameMap, infix(inParam, s, e));
							CommandLineOption const & opt = me.optionMap[option_index];
							s = --e;
							if (isBooleanOption(opt))
								_assignOptionValue(me, option_index, "true", estream);
							else
							{
								if (i + opt.argumentsPerOption < argc)
								{
									for (int t = 0; t < opt.argumentsPerOption; ++t)
										if (!_assignOptionValue(me, option_index, argv[++i], t, estream)) return false;
								}
								else // no value available
								{
									_streamWrite(estream, me.appName);
									_streamWrite(estream, ": \'");
									_writeOptName(estream, opt);
									_streamWrite(estream, "\' requires ");
									_streamPutInt(estream, opt.argumentsPerOption);
									_streamWrite(estream, " value(s)\n");
									return false;
								}
							}
						}
					}
					if (s == e)
					{
                        _streamWrite(estream, me.appName);
                        _streamWrite(estream, ": invalid option '-");
                        _streamWrite(estream, suffix(inParam, s));
                        _streamWrite(estream, "\'\n");
                        return false;
					}
				}
            }
            else if (inParam[1] == '-') // this is a long option
            {
                unsigned t = 2;
                CharString longOpt, val;
                for (; t < len && inParam[t] != '='; ++t)
					appendValue(longOpt, inParam[t], Generous());
                if (t < len) // this one is a --name=value option
					val = suffix(inParam, t + 1);
				
                // we may be got already a value
                if (hasOptionLong(me, longOpt))
                {
                    TOptionPosition option_index = cargo(me.longNameMap, longOpt);
                    CommandLineOption opt = me.optionMap[option_index];

                    if (!empty(val))
                    {
						if (opt.argumentsPerOption == 1)
						{
							if (!_assignOptionValue(me, option_index, val, estream)) return false;
						}
						else
						{
							_streamWrite(estream, me.appName);
							_streamWrite(estream, ": \'");
							_writeOptName(estream, opt);
							_streamWrite(estream, "\' requires ");
							_streamPutInt(estream, opt.argumentsPerOption);
							_streamWrite(estream, " values\n");
							return false;
						}
                    }
                    else if(isBooleanOption(opt))
                    {
						_assignOptionValue(me, option_index, "true", estream);
                    }
                    else if (i + opt.argumentsPerOption < argc)
					{
						for (int t = 0; t < opt.argumentsPerOption; ++t)
							if (!_assignOptionValue(me, option_index, argv[++i], t, estream)) return false;
					}
					else // no value available
					{
						_streamWrite(estream, me.appName);
						_streamWrite(estream, ": \'");
						_writeOptName(estream, opt);
						_streamWrite(estream, "\' requires ");
						_streamPutInt(estream, opt.argumentsPerOption);
						_streamWrite(estream, " value(s)\n");
						return false;
					}
                }
                else
                {
                    _streamWrite(estream, me.appName);
                    _streamWrite(estream, ": invalid option \'--");
                    _streamWrite(estream, longOpt);
                    _streamWrite(estream, "'\n");
                    return false;
                }
            }            
        }
        else
        { // this seems to be a normal argument
            appendValue(me.arguments,argv[i] );
        }
    }
	if (isSetLong(me, "version"))
	{
		version(me, estream);
        return true;
	}
    if (isSetLong(me, "help"))
    {
        help(me, estream);
        return true;
    }
	return _allMandatorySet(me) && (length(me.arguments) >= me.required_arguments);
}

//////////////////////////////////////////////////////////////////////////////

inline bool
parse(CommandLineParser & me, int argc, const char *argv[])
{
    return parse(me, argc, argv, ::std::cerr);
}


//////////////////////////////////////////////////////////////////////////////

inline String<CharString> const &
getOptionValues(CommandLineParser const & me, unsigned option_index)
{
    CommandLineOption const & opt = me.optionMap[option_index];
	if (empty(opt.value))
		return opt.defaultValue;
	else
		return opt.value;
}

inline CharString const &
getOptionValue(CommandLineParser const & me, unsigned option_index, unsigned argNo)
{
    CommandLineOption const & opt = me.optionMap[option_index];
	if (argNo < length(opt.value))
		return opt.value[argNo];
	else if (argNo < length(opt.defaultValue))
		return opt.defaultValue[argNo];
	else
		return me.null;
}

inline CharString const &
getOptionValue(CommandLineParser const & me, unsigned option_index)
{
	return getOptionValue(me, option_index, 0);
}

//////////////////////////////////////////////////////////////////////////////

inline bool
_convertOptionValue(CommandLineOption const & opt, bool & dst, CharString const & src)
{
    if (!isBooleanOption(opt)) return false;
	dst = !empty(src);
	return true;
}

inline bool
_convertOptionValue(CommandLineOption const & opt, int & dst, CharString const & src)
{
    if (!isIntOption(opt)) return false;
	dst = atoi(toCString(src));
	return length(src) > 0;
}

inline bool
_convertOptionValue(CommandLineOption const & opt, unsigned int & dst, CharString const & src)
{
    if (!isIntOption(opt)) return false;
	dst = atoi(toCString(src));
	return length(src) > 0;
}

inline bool
_convertOptionValue(CommandLineOption const & opt, float & dst, CharString const & src)
{
    if (!isDoubleOption(opt)) return false;
	dst = (float)atof(toCString(src));
	return length(src) > 0;
}

inline bool
_convertOptionValue(CommandLineOption const & opt, double & dst, CharString const & src)
{
    if (!isDoubleOption(opt)) return false;

	dst = atof(toCString(src));
	return length(src) > 0;
}

template <typename TObject>
inline bool
_convertOptionValue(CommandLineOption const & opt, TObject & dst, CharString const & src)
{
    if (!isStringOption(opt)) return false;
	assign(dst, src);
	return true;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getOptionValueShort:
..summary:Fills the passed variable $value$ with the value set for the option on the command line.
..cat:Miscellaneous
..signature:getOptionValueShort(parser,optionIdentifier,value)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the short option.
..param.value:The variable where the value is stored.
...remarks: The variable type ($int$, $double$, $bool$ or @Shortcut.CharString@) depends on the OptionTyp.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
*/
template <typename TValue>
inline bool
getOptionValueShort(CommandLineParser & me, CharString const & shortName, unsigned argNo, TValue & val)
{
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if (!hasOptionShort(me, shortName))
	{
		_streamWrite(std::cerr, me.appName);
		_streamWrite(std::cerr, ": \'");
		_streamWrite(std::cerr, shortName);
		_streamWrite(std::cerr, "\' is not an option\n");
		return false;
	}
    TOptionPosition option_index = cargo(me.shortNameMap, shortName);
	return _convertOptionValue(me.optionMap[option_index], val, getOptionValue(me, option_index, argNo));
}

template <typename TValue>
inline bool
getOptionValueShort(CommandLineParser & me, CharString const & shortName, TValue & val)
{
	return getOptionValueShort(me, shortName, 0, val);
}

inline String<CharString> const &
getOptionValuesShort(CommandLineParser & me,CharString const & shortName)
{
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if (!hasOptionShort(me, shortName))
	{
		_streamWrite(std::cerr, me.appName);
		_streamWrite(std::cerr, ": \'");
		_streamWrite(std::cerr, shortName);
		_streamWrite(std::cerr, "\' is not an option\n");
		return me.nullSet;
	}
    TOptionPosition option_index = cargo(me.shortNameMap, shortName);
	return getOptionValues(me, option_index);
}


/**
.Function.getOptionValueLong:
..summary:Fills the passed variable $value$ with the value set for the option on the command line.
..cat:Miscellaneous
..signature:getOptionValueLong(parser,optionIdentifier,value)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the short option.
..param.value:The variable where the value is stored.
...remarks: The variable type ($int$, $double$, $bool$ or @Shortcut.CharString@) depends on the OptionTyp.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
*/
template <typename TValue>
inline bool
getOptionValueLong(CommandLineParser & me, CharString const & longName, unsigned argNo, TValue & val)
{
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if (!hasOptionLong(me, longName))
	{
		_streamWrite(std::cerr, me.appName);
		_streamWrite(std::cerr, ": \'");
		_streamWrite(std::cerr, longName);
		_streamWrite(std::cerr, "\' is not an option\n");
		return false;
	}
    TOptionPosition option_index = cargo(me.longNameMap, longName);
	return _convertOptionValue(me.optionMap[option_index], val, getOptionValue(me, option_index, argNo));
}

template <typename TValue>
inline bool
getOptionValueLong(CommandLineParser & me,CharString const & longName, TValue & val)
{
	return getOptionValueLong(me, longName, 0, val);
}

inline String<CharString> const &
getOptionValuesLong(CommandLineParser & me, CharString const & longName)
{
    typedef Size<String<CommandLineOption> >::Type TOptionPosition;

    if (!hasOptionLong(me, longName))
	{
		_streamWrite(std::cerr, me.appName);
		_streamWrite(std::cerr, ": \'");
		_streamWrite(std::cerr, longName);
		_streamWrite(std::cerr, "\' is not an option\n");
		return me.nullSet;
	}
    TOptionPosition option_index = cargo(me.longNameMap, longName);
	return getOptionValues(me, option_index);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getArgumentValue:
..summary:Fills the passed variable $value$ with the argument set on the command line.
..cat:Miscellaneous
..signature:getArgumentValue(parser,position,value)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.position:A zero based $int$ indicating which argument you want to get.
..param.value:The variable where the value is stored.
...type:Shortcut.CharString
..returns: $true$ if the requested argument exists, $false$ otherwise.
*/
inline CharString const &
getArgumentValue(CommandLineParser const & me, unsigned position)
{
    if (position < length(me.arguments))
        return me.arguments[position];
    else
		return me.null;
}

inline String<CharString> const &
getArgumentValues(CommandLineParser const & me)
{
	return me.arguments;
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.argumentCount:
..summary:Returns the count of passed arguments.
..cat:Miscellaneous
..signature:argumentCount(parser)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
*/
inline Size<String<CharString> >::Type
argumentCount(CommandLineParser const & me)
{
    return length(me.arguments);
}


} // end SEQAN_NAMESPACE_MAIN

#endif
