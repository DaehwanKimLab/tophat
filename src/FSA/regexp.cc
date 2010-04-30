// Regexp.cc - regular expression class based on Henry Spencer's regexp code.
// Adapted for Windows by Guy Gascoigne (see notes below)
// Made unix-friendly by Ian Holmes ihh@lanl.gov April 30, 1999

// In case this isn't obvious from the later comments this is an ALTERED
// version of the software. If you like my changes then cool, but nearly
// all of the functionality here is derived from Henry Spencer's original
// work.
//
// This code should work correctly under both _SBCS and _UNICODE, I did
// start working on making it work with _MBCS but gave up after a while
// since I don't need this particular port and it's not going to be as
// straight forward as the other two.
//
// The problem stems from the compiled program being stored as chars,
// the individual items need to be wide enough to hold whatever character
// is thrown at them, but currently they are accessed as an array of
// whatever size integral type is appropriate.  _MBCS would cause this
// to be char, but at times it would need to be larger.  This would
// require making the program be an array of short with the appropriate
// conversions used everywhere.  Certainly it's doable, but it's a pain.
// What's worse is that the current code will compile and run under _MBCS,
// only breaking when it gets wide characters thrown against it.
//
// I've marked at least one bit of code with #pragma messages, I may not
// get all of them, but they should be a start
//
// Guy Gascoigne - Piggford (ggp@bigfoot.com) Friday, February 27, 1998

// EAT MY SHORTS - I took out all the casts to shorts cos they made my
// SparcStation unhappy. ihh@lanl.gov April 30, 1999

// regcomp and regexec -- regsub and regerror are elsewhere
// @(#)regexp.c	1.3 of 18 April 87
//
//	Copyright (c) 1986 by University of Toronto.
//	Written by Henry Spencer.  Not derived from licensed software.
//
//	Permission is granted to anyone to use this software for any
//	purpose on any computer system, and to redistribute it freely,
//	subject to the following restrictions:
//
//	1. The author is not responsible for the consequences of use of
//		this software, no matter how awful, even if they arise
//		from defects in it.
//
//	2. The origin of this software must not be misrepresented, either
//		by explicit claim or by omission.
//
//	3. Altered versions must be plainly marked as such, and must not
//		be misrepresented as being the original software.
// *** THIS IS AN ALTERED VERSION.  It was altered by John Gilmore,
// *** hoptoad!gnu, on 27 Dec 1986, to add \< and \> for word-matching
// *** as in BSD grep and ex.
// *** THIS IS AN ALTERED VERSION.  It was altered by John Gilmore,
// *** hoptoad!gnu, on 28 Dec 1986, to optimize characters quoted with \.
// *** THIS IS AN ALTERED VERSION.  It was altered by James A. Woods,
// *** ames!jaw, on 19 June 1987, to quash a regcomp() redundancy.
// *** THIS IS AN ALTERED VERSION.  It was altered by Geoffrey Noer,
// *** THIS IS AN ALTERED VERSION.  It was altered by Guy Gascoigne - Piggford
// *** guy@wyrdrune.com, on 15 March 1998, porting it to C++ and converting
// *** it to be the engine for the Regexp class
// *** THIS IS AN ALTERED VERSION.  It was altered by Ian Holmes,
// *** ihh@fruitfly.org, on Jan 28 2000, porting it to gcc and DART's string
// *** class (sstring)
// *** THIS IS AN ALTERED VERSION.  It was altered by Robert Bradley,
// *** Dec 4 2008, changing the string class: sstring -> std::string.
//
// Beware that some of this code is subtly aware of the way operator
// precedence is structured in regular expressions.  Serious changes in
// regular-expression syntax might require a total rethink.

#include <cassert>
#include <iostream>
#include <cstring>

#include "regexp.h"

// The first byte of the regexp internal "program" is actually this magic
// number; the start node begins in the second byte.

const char	MAGIC = ((char)'\234');

// #define new DEBUG_NEW

// following line commented out by ihh, December 14 2003, because gcc doesn't care for it
// #pragma warning( disable : 4711 )	// automatic inline selected

// The "internal use only" fields in regexp.h are present to pass info from
// compile to execute that permits the execute phase to run lots faster on
// simple cases.  They are:
//
// regstart	char that must begin a match; '\0' if none obvious
// reganch	is the match anchored (at beginning-of-line only)?
// regmust	string (pointer into program) that match must include, or NULL
// regmlen	length of regmust string
//
// Regstart and reganch permit very fast decisions on suitable starting
// points for a match, cutting down the work a lot.  Regmust permits fast
// rejection of lines that cannot possibly match.  The regmust tests are
// costly enough that regcomp() supplies a regmust only if the
// r.e. contains something potentially expensive (at present, the only
// such thing detected is * or + at the start of the r.e., which can
// involve a lot of backup).  Regmlen is supplied because the test in
// regexec() needs it and regcomp() is computing it anyway.

// Structure for regexp "program".  This is essentially a linear encoding
// of a nondeterministic finite-state machine (aka syntax charts or
// "railroad normal form" in parsing technology).  Each node is an opcode
// plus a "next" pointer, possibly plus an operand.  "Next" pointers of
// all nodes except BRANCH implement concatenation; a "next" pointer with
// a BRANCH on both ends of it is connecting two alternatives.  (Here we
// have one of the subtle syntax dependencies: an individual BRANCH (as
// opposed to a collection of them) is never concatenated with anything
// because of operator precedence.)  The operand of some types of node is
// a literal string; for others, it is a node leading into a sub-FSM.  In
// particular, the operand of a BRANCH node is the first node of the
// branch.  (NB this is *not* a tree structure: the tail of the branch
// connects to the thing following the set of BRANCHes.)  The opcodes
// are:

enum {
//  definition	number		opnd?	meaning
	END		=	0,		//	no		End of program.
	BOL		=	1,		//	no		Match beginning of line.
	EOL		=	2,		//	no		Match end of line.
	ANY		=	3,		//	no		Match any character.
	ANYOF	=	4,		//	str		Match any of these.
	ANYBUT	=	5,		//	str		Match any but one of these.
	BRANCH	=	6,		//	node	Match this, or the next..\&.
	BACK	=	7,		//	no		"next" ptr points backward.
	EXACTLY	=	8,		//	str		Match this string.
	NOTHING	=	9,		//	no		Match empty string.
	STAR	=	10,		//	node	Match this 0 or more times.
	PLUS	=	11,		//	node	Match this 1 or more times.
	WORDA	=	12,		//	no		Match "" at wordchar, where prev is nonword
	WORDZ	=	13,		//	no		Match "" at nonwordchar, where prev is word
	OPEN	=	20,		//	no		Sub-RE starts here.
						//			OPEN+1 is number 1, etc.
	CLOSE	=	30		//	no		Analogous to OPEN.
};

// Opcode notes:
//
// BRANCH	The set of branches constituting a single choice are hooked
//		together with their "next" pointers, since precedence prevents
//		anything being concatenated to any individual branch.  The
//		"next" pointer of the last BRANCH in a choice points to the
//		thing following the whole choice.  This is also where the
//		final "next" pointer of each individual branch points; each
//		branch starts with the operand node of a BRANCH node.
//
// BACK		Normal "next" pointers all implicitly point forward; BACK
//		exists to make loop structures possible.
//
// STAR,PLUS	'?', and complex '*' and '+', are implemented as circular
//		BRANCH structures using BACK.  Simple cases (one character
//		per match) are implemented with STAR and PLUS for speed
//		and to minimize recursive plunges.
//
// OPEN,CLOSE	...are numbered at compile time.

// A node is one char of opcode followed by two chars of "next" pointer.
// "Next" pointers are stored as two 8-bit pieces, high order first.  The
// value is a positive offset from the opcode of the node containing it.
// An operand, if any, simply follows the node.  (Note that much of the
// code generation knows about this implicit relationship.)
//
// Using two bytes for the "next" pointer is vast overkill for most things,
// but allows patterns to get big without disasters.


enum
{
	REGERR_SENTINEL_VALUE = 0,
	REGERR_NULLARG = 1, REGERR_CORRUPTED, REGERR_CORRUPTION, REGERR_CORRUPTED_POINTERS,
	REGERR_BAD_REGREPEAT, REGERR_CORRUPTED_OPCODE, REGERR_NULL_TO_REGSUB,
	REGERR_DAMAGED_REGEXP_REGSUB, REGERR_DAMAGED_MATCH_STRING, REGERR_NULL_TO_REGCOMP,
	REGERR_TO_BIG, REGERR_TO_MANY_PAREN, REGERR_UNTERMINATED_PAREN, REGERR_UNMATCHED_PAREN,
	REGERR_INTERNAL_ERROR_JUNK, REGERR_OP_COULD_BE_EMPTY, REGERR_NESTED_OP, REGERR_INVALID_RANGE,
	REGERR_UNMATCHED_BRACE, REGERR_INTERNAL_UNEXPECTED_CHAR, REGERR_OP_FOLLOWS_NOTHING,
	REGERR_TRAILING_ESC, REGERR_INTERNAL_STRSCSPN, REGERR_NO_REGEXP
};

struct regErr
{
	int m_id;
	const char* m_err;
} errors[] = {
	{ REGERR_NULLARG,					"NULL argument to regexec" },
	{ REGERR_CORRUPTED,					"corrupted regexp" },
	{ REGERR_CORRUPTION,				"regexp corruption" },
	{ REGERR_CORRUPTED_POINTERS,		"corrupted pointers" },
	{ REGERR_BAD_REGREPEAT,				"internal error: bad call of regrepeat" },
	{ REGERR_CORRUPTED_OPCODE,			"corrupted opcode" },
	{ REGERR_NULL_TO_REGSUB,			"NULL parm to regsub" },
	{ REGERR_DAMAGED_REGEXP_REGSUB,		"damaged regexp fed to regsub" },
	{ REGERR_DAMAGED_MATCH_STRING,		"damaged match string" },
	{ REGERR_NULL_TO_REGCOMP,			"NULL argument to regcomp" },
	{ REGERR_TO_BIG,					"regexp too big" },
	{ REGERR_TO_MANY_PAREN,				"too many ()" },
	{ REGERR_UNTERMINATED_PAREN,		"unterminated ()" },
	{ REGERR_UNMATCHED_PAREN,			"unmatched ()" },
	{ REGERR_INTERNAL_ERROR_JUNK,		"internal error: junk on end" },
	{ REGERR_OP_COULD_BE_EMPTY,			"*+ operand could be empty" },
	{ REGERR_NESTED_OP,					"nested *?+" },
	{ REGERR_INVALID_RANGE,				"invalid [] range" },
	{ REGERR_UNMATCHED_BRACE,			"unmatched []" },
	{ REGERR_INTERNAL_UNEXPECTED_CHAR,	"internal error: \\0|) unexpected" },
	{ REGERR_OP_FOLLOWS_NOTHING,		"?+* follows nothing" },
	{ REGERR_TRAILING_ESC,				"trailing \\" },
	{ REGERR_INTERNAL_STRSCSPN,			"internal error: strcspn 0" },
	{ REGERR_NO_REGEXP,					"NULL regexp" },
	{ REGERR_SENTINEL_VALUE,			"Unknown error" }							// must be last value
};

// Flags to be passed up and down.

enum {
	HASWIDTH	=	01,	// Known never to match null string.
	SIMPLE		=	02,	// Simple enough to be STAR/PLUS operand.
	SPSTART		=	04,	// Starts with * or +.
	WORST		=	0	// Worst case.
};

///////////////////////////////////////////////////////////////////////////////

class CRegErrorHandler
{
	friend class Regexp;
	mutable std::string m_szError;
	static const char* FindErr( int id );
protected:
	void ClearErrorString() const;
	void regerror( const char* s ) const;
	void regerror( int id ) const;
public:
	CRegErrorHandler() { }
	CRegErrorHandler( const CRegErrorHandler & reh ) : m_szError( reh.m_szError ) {}

	const std::string & GetErrorString() const;
};

void CRegErrorHandler::regerror( const char* s ) const
{
  std::cerr << "regerror: " << s << "\n";
  m_szError = s;
}

void CRegErrorHandler::regerror( int id ) const
{
	regerror( FindErr( id ) );
}

const std::string & CRegErrorHandler::GetErrorString() const
{
	return m_szError;
}

void CRegErrorHandler::ClearErrorString() const
{
	m_szError.clear();
}

const char* CRegErrorHandler::FindErr( int id )
{
  struct regErr * perr;
  for (perr = errors; perr->m_id != REGERR_SENTINEL_VALUE; perr++ )
    if ( perr->m_id == id )
      return perr->m_err;
  
  return perr->m_err;		// since we've fallen off the array, perr->m_id == 0
}

///////////////////////////////////////////////////////////////////////////////

// All of the functions required to directly access the 'program'
class CRegProgramAccessor : public CRegErrorHandler
{
public:
	static inline char OP( char* p )
	{
		return (*(p));
	}
	static inline char* OPERAND( char* p )
	{
		return p + 3;
	}
	static inline char* regnext( char* p )
	{
	  const int offset = (((unsigned char*)p)[1] << 8) + ((unsigned char*)p)[2];
	  
		if (offset == 0)
			return(NULL);
		
		return((OP(p) == BACK) ? p-offset : p+offset);
	}
#ifdef _RE_DEBUG
	char* CRegProgramAccessor::regprop( char* op );
#endif
};

///////////////////////////////////////////////////////////////////////////////

// The internal interface to the regexp, wrapping the compilation as well as the
// execution of the regexp (matching)

class regexp_internal : public CRegProgramAccessor
{
	friend class CRegExecutor;
	friend class Regexp;
	
	int m_programSize;
	char* startp[Regexp::NSUBEXP + 1];
	char* endp[Regexp::NSUBEXP + 1];
	char regstart;		// Internal use only.
	char reganch;		// Internal use only.
	char* regmust;		// Internal use only.
	int regmlen;		// Internal use only.
	char* program;
	
	bool status;
	int count;			// used by Regexp to manage the reference counting of regexps
	int numSubs;
public:
	
	regexp_internal( const char* exp, bool iCase );
	regexp_internal( const regexp_internal & r );
	~regexp_internal();
	
	void ignoreCase( const char * in, char * out );
	
	bool regcomp( const char* exp );
	bool regexec( char* str );
	bool Status() const { return status; }

	std::string GetReplaceString( const char* sReplaceExp ) const;

	regexp_internal * getCopy();

#ifdef _RE_DEBUG
	void  regdump();
#endif
	
#ifdef _DEBUG
	std::string m_originalPattern;
	std::string m_modifiedPattern;
#endif
};

///////////////////////////////////////////////////////////////////////////////
// Compile / Validate the regular expression - ADT

class CRegCompilerBase : public CRegProgramAccessor
{
public:
	CRegCompilerBase( const char* parse );

	char* reg(int paren, int *flagp);
protected:
	char* regparse;	// Input-scan pointer. 
	int regnpar;		// () count. 

	char* regbranch(int *flagp);
	char* regpiece(int *flagp);
	char* regatom(int *flagp);
	inline bool ISREPN( char c )	{ return ((c) == '*' || (c) == '+' || (c) == '?'); }

	virtual void regc(int c) = 0;
	virtual char* regnode(int op) = 0;
	virtual void reginsert(char op, char* opnd) = 0;
	virtual void regtail(char* p, char* val) = 0;
	virtual void regoptail(char* p, char* val) = 0;

  // virtual destructor, to keep gcc happy -- added by ihh, December 14, 2003
  virtual ~CRegCompilerBase() { }
};

///////////////////////////////////////////////////////////////////////////////
// First pass over the expression, testing for validity and returning the 
// program size

class CRegValidator : public CRegCompilerBase
{
public:
	CRegValidator( const char* parse );

	inline long Size() const { return regsize; }
private:
	long regsize;		// Code size. 
	char regdummy[3];	// NOTHING, 0 next ptr 
protected:
	char* regnode(int) { regsize += 3; return regdummy; }
	void regc(int) { regsize++; }
	void reginsert(char, char*) { regsize += 3; }
	void regtail(char*, char*) { return; }
	void regoptail(char*, char*) { return; }
};

///////////////////////////////////////////////////////////////////////////////
// Second pass, actually generating the program

class CRegCompiler : public CRegCompilerBase
{
public:
	CRegCompiler( const char* parse, char* prog );
private:
	char* regcode;
protected:
	// regc - emit (if appropriate) a byte of code
	void regc(int b)
	{
		*regcode++ = (char)b;
	}
	char* regnode(int op);
	void reginsert(char op, char* opnd);
	void regtail(char* p, char* val);
	void regoptail(char* p, char* val);
};

// regnode - emit a node
char* CRegCompiler::regnode(int op)
{
	char* const ret = regcode;

	char* ptr = ret;
	*ptr++ = (char)op;
	*ptr++ = '\0';		// Null next pointer. 
	*ptr++ = '\0';
	regcode = ptr;

	return(ret);
}

// reginsert - insert an operator in front of already-emitted operand
//
// Means relocating the operand.
void CRegCompiler::reginsert(char op, char* opnd)
{
	char* place;

	(void) memmove(opnd+3, opnd, (size_t)((regcode - opnd)*sizeof(char)));
	regcode += 3;

	place = opnd;		// Op node, where operand used to be.
	*place++ = op;
	*place++ = '\0';
	*place++ = '\0';
}

// regtail - set the next-pointer at the end of a node chain
void CRegCompiler::regtail(char* p, char* val)
{
	char* scan;
	char* temp;

	// Find last node. 
	for (scan = p; (temp = regnext(scan)) != NULL; scan = temp)
		continue;

	int next_ptr = (OP(scan) == BACK) ? scan - val : val - scan;
	((unsigned char*)scan)[1] = (unsigned char) (next_ptr >> 8);
	((unsigned char*)scan)[2] = (unsigned char) (next_ptr & 0xff);

	// debugging assertion added by ihh -- December 5, 2003
	if (regnext(scan) != val)
	  regerror( REGERR_CORRUPTED );
}

// regoptail - regtail on operand of first argument; nop if operandless
void CRegCompiler::regoptail(char* p, char* val)
{
	// "Operandless" and "op != BRANCH" are synonymous in practice. 
	if (OP(p) == BRANCH)
		regtail(OPERAND(p), val);
}

///////////////////////////////////////////////////////////////////////////////

CRegCompilerBase::CRegCompilerBase( const char* parse )
	: regparse( (char*)parse ),
	regnpar(1)
{
}

CRegValidator::CRegValidator( const char* parse )
	: CRegCompilerBase( parse ),
	regsize(0)
{
	regc(MAGIC);
	regdummy[0] = NOTHING;
	regdummy[1] = regdummy[2] = 0;
}

CRegCompiler::CRegCompiler( const char* parse, char* prog )
	: CRegCompilerBase( parse ),
	regcode(prog)
{
	regc(MAGIC);
}

///////////////////////////////////////////////////////////////////////////////

regexp_internal::regexp_internal( const char* exp, bool iCase )
  : m_programSize(0),
    regstart(0),
    reganch(0),
    regmust(0),
    regmlen(0),
    program(0)
{
#if _DEBUG
	m_originalPattern = exp;		// keep a version of the pattern for debugging
#endif

	if ( iCase )
	{
		char * out = new char[(strlen( exp ) * 4) + 1];
		ignoreCase( exp, out );

#if _DEBUG
		m_modifiedPattern = out;	// and the modified version if there is one
#endif
		status = regcomp( out );
		delete [] out;
	}
	else
		status = regcomp( exp );

	count = numSubs = 0;
}

regexp_internal::regexp_internal( const regexp_internal & orig )
	: m_programSize(orig.m_programSize),
	  regstart(orig.regstart),
	  reganch(orig.reganch),
	  regmust(0),
	  regmlen(orig.regmlen),
	  numSubs(orig.numSubs)
{
#if _DEBUG
	m_originalPattern = orig.m_originalPattern;
	m_modifiedPattern = orig.m_modifiedPattern;
#endif
	status = orig.status;
	count = 0;
	program = new char[m_programSize];
	memcpy( program, orig.program, m_programSize * sizeof( char ) );
	if ( orig.regmust )
		regmust = program + ( orig.regmust - orig.program );

	for ( int i = Regexp::NSUBEXP; i > 0; i--)
	{
		startp[i] = orig.startp[i];
		endp[i] = orig.endp[i];
	}
}

regexp_internal::~regexp_internal()
{
	delete [] program;
}


// regcomp - compile a regular expression into internal code
//
// We can't allocate space until we know how big the compiled form will
// be, but we can't compile it (and thus know how big it is) until we've
// got a place to put the code.  So we cheat: we compile it twice, once
// with code generation turned off and size counting turned on, and once
// "for real".  This also means that we don't allocate space until we are
// sure that the thing really will compile successfully, and we never
// have to move the code and thus invalidate pointers into it.  (Note
// that it has to be in one piece because free() must be able to free it
// all.)
//
// Beware that the optimization-preparation code in here knows about some
// of the structure of the compiled regexp.

bool regexp_internal::regcomp(const char* exp)
{
	char* scan;
	int flags;

	if (exp == NULL)
	{
		regerror( REGERR_NULL_TO_REGCOMP );
		return 0;
	}

	// First pass: determine size, legality.
	CRegValidator tester( exp );

	if (tester.reg(0, &flags) == NULL)
		return false;

	// Small enough for pointer-storage convention? 
	if (tester.Size() >= 0x7fffL)	// Probably could be 0xffffL. 
	{
		regerror(REGERR_TO_BIG);
		return 0;
	}

	m_programSize = tester.Size();
	// Allocate space. 
	program = new char[m_programSize];

	CRegCompiler comp( exp, program );
	// Second pass: emit code. 
	if (comp.reg(0, &flags) == NULL)
		return false;

	scan = program + 1;			// First BRANCH. 
	if (OP(regnext(scan)) == END)
	{		// Only one top-level choice. 
		scan = OPERAND(scan);

		// Starting-point info. 
		if (OP(scan) == EXACTLY)
			regstart = *OPERAND(scan);
		else if (OP(scan) == BOL)
			reganch = 1;

		// If there's something expensive in the r.e., find the
		// longest literal string that must appear and make it the
		// regmust.  Resolve ties in favor of later strings, since
		// the regstart check works with the beginning of the r.e.
		// and avoiding duplication strengthens checking.  Not a
		// strong reason, but sufficient in the absence of others.
		 
		if (flags&SPSTART) 
		{
			char* longest = NULL;
			size_t len = 0;

			for (; scan != NULL; scan = regnext(scan))
				if (OP(scan) == EXACTLY && strlen(OPERAND(scan)) >= len)
				{
					longest = OPERAND(scan);
					len = strlen(OPERAND(scan));
				}
			regmust = longest;
			regmlen = (int)len;
		}
	}

	return true;
}

regexp_internal * regexp_internal::getCopy()
{
	return new regexp_internal( *this );
}

// reg - regular expression, i.e. main body or parenthesized thing
//
// Caller must absorb opening parenthesis.
//
// Combining parenthesis handling with the base level of regular expression
// is a trifle forced, but the need to tie the tails of the branches to what
// follows makes it hard to avoid.

char* CRegCompilerBase::reg( int paren, int *flagp )
{
	char* ret = 0;
	char* br;
	char* ender;
	int parno = 0;
	int flags;

	*flagp = HASWIDTH;	// Tentatively. 

	if (paren)
	{
		// Make an OPEN node. 
		if (regnpar >= Regexp::NSUBEXP)
		{
			regerror(REGERR_TO_MANY_PAREN);
			return NULL;
		}
		parno = regnpar;
		regnpar++;
		ret = regnode(OPEN+parno);
	}

	// Pick up the branches, linking them together. 
	br = regbranch(&flags);
	if (br == NULL)
		return(NULL);
	if (paren)
		regtail(ret, br);	// OPEN -> first. 
	else
		ret = br;
	*flagp &= ~(~flags&HASWIDTH);	// Clear bit if bit 0. 
	*flagp |= flags&SPSTART;
	while (*regparse == '|')
	{
		regparse++;
		br = regbranch(&flags);
		if (br == NULL)
			return(NULL);
		regtail(ret, br);	// BRANCH -> BRANCH. 
		*flagp &= ~(~flags&HASWIDTH);
		*flagp |= flags&SPSTART;
	}

	// Make a closing node, and hook it on the end. 
	ender = regnode((paren) ? CLOSE+parno : END);
	regtail(ret, ender);

	// Hook the tails of the branches to the closing node. 
	for (br = ret; br != NULL; br = regnext(br))
		regoptail(br, ender);

	// Check for proper termination. 
	if (paren && *regparse++ != ')')
	{
		regerror( REGERR_UNTERMINATED_PAREN );
		return NULL;
	}
	else if (!paren && *regparse != '\0')
	{
		if (*regparse == ')')
		{
			regerror( REGERR_UNMATCHED_PAREN );
			return NULL;
		}
		else
		{
			regerror( REGERR_INTERNAL_ERROR_JUNK );
			return NULL;
		}
		// NOTREACHED 
	}

	return(ret);
}

// regbranch - one alternative of an | operator
//
// Implements the concatenation operator.

char* CRegCompilerBase::regbranch(int *flagp)
{
	char* ret;
	char* chain;
	char* latest;
	int flags;
	int c;

	*flagp = WORST;				// Tentatively. 

	ret = regnode(BRANCH);
	chain = NULL;
	while ((c = *regparse) != '\0' && c != '|' && c != ')')
	{
		latest = regpiece(&flags);
		if (latest == NULL)
			return(NULL);
		*flagp |= flags&HASWIDTH;
		if (chain == NULL)		// First piece. 
			*flagp |= flags&SPSTART;
		else
			regtail(chain, latest);
		chain = latest;
	}
	if (chain == NULL)			// Loop ran zero times. 
		(void) regnode(NOTHING);

	return(ret);
}

// regpiece - something followed by possible [*+?]
//
// Note that the branching code sequences used for ? and the general cases
// of * and + are somewhat optimized:  they use the same NOTHING node as
// both the endmarker for their branch list and the body of the last branch.
// It might seem that this node could be dispensed with entirely, but the
// endmarker role is not redundant.

char* CRegCompilerBase::regpiece(int *flagp)
{
	char* ret;
	char op;
	char* next;
	int flags;

	ret = regatom(&flags);
	if (ret == NULL)
		return(NULL);

	op = *regparse;
	if (!ISREPN(op))
	{
		*flagp = flags;
		return(ret);
	}

	if (!(flags&HASWIDTH) && op != '?' )
	{
		regerror( REGERR_OP_COULD_BE_EMPTY );
		return NULL;
	}
	switch (op)
	{
		case '*':	*flagp = WORST|SPSTART;				break;
		case '+':	*flagp = WORST|SPSTART|HASWIDTH;	break;
		case '?':	*flagp = WORST;						break;
	}

	if (op == '*' && (flags&SIMPLE))
		reginsert(STAR, ret);
	else if (op == '*')
	{
		// Emit x* as (x&|), where & means "self". 
		reginsert(BRANCH, ret);				// Either x 
		regoptail(ret, regnode(BACK));		// and loop 
		regoptail(ret, ret);				// back 
		regtail(ret, regnode(BRANCH));		// or 
		regtail(ret, regnode(NOTHING));		// null. 
	}
	else if (op == '+' && (flags&SIMPLE))
		reginsert(PLUS, ret);
	else if (op == '+' )
	{
		// Emit x+ as x(&|), where & means "self". 
		next = regnode(BRANCH);				// Either 
		regtail(ret, next);
		regtail(regnode(BACK), ret);		// loop back 
		regtail(next, regnode(BRANCH));		// or 
		regtail(ret, regnode(NOTHING));		// null. 
	}
	else if (op == '?' )
	{
		// Emit x? as (x|) 
		reginsert(BRANCH, ret);				// Either x 
		regtail(ret, regnode(BRANCH));		// or 
		next = regnode(NOTHING);			// null. 
		regtail(ret, next);
		regoptail(ret, next);
	}
	regparse++;
	if (ISREPN(*regparse))
	{
		regerror( REGERR_NESTED_OP );
		return NULL;
	}

	return(ret);
}

// regatom - the lowest level
//
// Optimization:  gobbles an entire sequence of ordinary characters so that
// it can turn them into a single node, which is smaller to store and
// faster to run.  Backslashed characters are exceptions, each becoming a
// separate node; the code is simpler that way and it's not worth fixing.

char* CRegCompilerBase::regatom(int * flagp)
{
	char* ret;
	int flags;

	*flagp = WORST;		// Tentatively. 

	switch ( *regparse++ )
	{
		// FIXME: these chars only have meaning at beg/end of pat?
		case '^':
			ret = regnode(BOL);
			break;
		case '$':
			ret = regnode(EOL);
			break;
		case '.':
			ret = regnode(ANY);
			*flagp |= HASWIDTH|SIMPLE;
			break;
		case '[':
		{
			int range;
			int rangeend;
			int c;

			if (*regparse == '^')
			{	// Complement of range. 
				ret = regnode(ANYBUT);
				regparse++;
			}
			else
				ret = regnode(ANYOF);
			if ((c = *regparse) == ']' || c == '-')
			{
				regc(c);
				regparse++;
			}
			while ((c = *regparse++ ) != '\0' && c != ']')
			{
				if (c != '-')
					regc(c);
				else if ((c = *regparse) == ']' || c == '\0')
					regc('-');
				else
				{
					range = (char)*(regparse-2);
					rangeend = (char)c;
					if (range > rangeend)
					{
						regerror( REGERR_INVALID_RANGE );
						return NULL;
					}
					for (range++; range <= rangeend; range++)
						regc(range);
					regparse++;
				}
			}
			regc('\0');
			if (c != ']')
			{
				regerror( REGERR_UNMATCHED_BRACE );
				return NULL;
			}
			*flagp |= HASWIDTH|SIMPLE;
			break;
		}
		case '(':
			ret = reg(1, &flags);
			if (ret == NULL)
				return(NULL);
			*flagp |= flags&(HASWIDTH|SPSTART);
			break;
		case '\0':
		case '|':
		case ')':
			// supposed to be caught earlier 
			regerror( REGERR_INTERNAL_UNEXPECTED_CHAR );
			return NULL;
		case '?':
		case '+':
		case '*':
			{
				regerror( REGERR_OP_FOLLOWS_NOTHING );
				return NULL;
			}
		case '\\':
			switch (*regparse++)
			{
				case '\0':
				{
					regerror( REGERR_TRAILING_ESC );
					return NULL;
				}
				case '<':
					ret = regnode(WORDA);
					break;
				case '>':
					ret = regnode(WORDZ);
					break;
				/* FIXME: Someday handle \1, \2, ... */
				default:
					/* Handle general quoted chars in exact-match routine */
					goto de_fault;
			}
			break;
	de_fault:
	default:
		// Encode a string of characters to be matched exactly.
		//
		// This is a bit tricky due to quoted chars and due to
		// '*', '+', and '?' taking the SINGLE char previous
		// as their operand.
		//
		// On entry, the char at regparse[-1] is going to go
		// into the string, no matter what it is.  (It could be
		// following a \ if we are entered from the '\' case.)
		// 
		// Basic idea is to pick up a good char in  ch  and
		// examine the next char.  If it's *+? then we twiddle.
		// If it's \ then we frozzle.  If it's other magic char
		// we push  ch  and terminate the string.  If none of the
		// above, we push  ch  on the string and go around again.
		//
		//  regprev  is used to remember where "the current char"
		// starts in the string, if due to a *+? we need to back
		// up and put the current char in a separate, 1-char, string.
		// When  regprev  is NULL,  ch  is the only char in the
		// string; this is used in *+? handling, and in setting
		// flags |= SIMPLE at the end.
		{
			char *regprev;
			register char ch;

			regparse--;			/* Look at cur char */
			ret = regnode(EXACTLY);
			for ( regprev = 0 ; ; ) {
				ch = *regparse++;	/* Get current char */
				switch (*regparse) {	/* look at next one */

				default:
					regc(ch);	/* Add cur to string */
					break;

				case '.': case '[': case '(':
				case ')': case '|': case '\n':
				case '$': case '^':
				case '\0':
				/* FIXME, $ and ^ should not always be magic */
				magic:
					regc(ch);	/* dump cur char */
					goto done;	/* and we are done */

				case '?': case '+': case '*':
					if (!regprev) 	/* If just ch in str, */
						goto magic;	/* use it */
					/* End mult-char string one early */
					regparse = regprev; /* Back up parse */
					goto done;

				case '\\':
					regc(ch);	/* Cur char OK */
					switch (regparse[1]){ /* Look after \ */
					case '\0':
					case '<':
					case '>':
					/* FIXME: Someday handle \1, \2, ... */
						goto done; /* Not quoted */
					default:
						/* Backup point is \, scan							 * point is after it. */
						regprev = regparse;
						regparse++; 
						continue;	/* NOT break; */
					}
				}
				regprev = regparse;	/* Set backup point */
			}
		done:
			regc('\0');
			*flagp |= HASWIDTH;
			if (!regprev)		/* One char? */
				*flagp |= SIMPLE;
		}
		break;
	}

	return(ret);
}

////////////////////////////////////////////////////////////////////////////////
// regexec and friends

// Work-variable struct for regexec().

class CRegExecutor : public CRegProgramAccessor
{
	friend bool regexp_internal::regexec( char* str );

	char* reginput;		// String-input pointer. 
	char* regbol;			// Beginning of input, for ^ check. 
	char* * regstartp;		// Pointer to startp array. 
	char* * regendp;		// Ditto for endp.

	regexp_internal * prog;
public:
	CRegExecutor( regexp_internal * prog, char* str );
protected:
	bool regtry( char* str );
	bool regmatch( char* prog );
	size_t regrepeat( char* node );
};

CRegExecutor::CRegExecutor( regexp_internal * p, char* str )
	: regbol( str ),
	regstartp( p->startp ),
	regendp( p->endp ),
	prog(p)
{
}

#ifdef _RE_DEBUG
int regnarrate = 0;
#endif

// regexec - match a regexp against a string

bool regexp_internal::regexec( char* str )
{
	// Be paranoid. 
	if ( str == NULL )
	{
		regerror( REGERR_NULLARG );
		return false;
	}

	// Check validity of program. 
	if (*program != MAGIC)
	{
		regerror( REGERR_CORRUPTED );
		return false;
	}

	// If there is a "must appear" string, look for it. 
	if ( regmust != NULL && strstr( str, regmust ) == NULL )
		return false;

	CRegExecutor executor( this, str );

	// Simplest case:  anchored match need be tried only once. 
	if ( reganch )
		return( executor.regtry( str ) );

	// Messy cases:  unanchored match. 
	if ( regstart != '\0' )
	{
		// We know what char it must start with. 
		for ( char* s = str; s != NULL; s = strchr( s+1 , regstart ) )
			if ( executor.regtry( s) )
				return true;
		return false;
	}
	else
	{
		// We don't -- general case. 
		for ( char* s = str; ! executor.regtry( s ); s++ )
			if (*s == '\0')
				return false;
	}
	return true;
}

// regtry - try match at specific point
bool CRegExecutor::regtry( char* str )
{
	int i;
	char* * stp;
	char* * enp;

	reginput = str;

	stp = prog->startp;
	enp = prog->endp;
	for (i = Regexp::NSUBEXP; i > 0; i--)
	{
		*stp++ = NULL;
		*enp++ = NULL;
	}
	if ( regmatch( prog->program + 1 ) )
	{
		prog->startp[0] = str;
		prog->endp[0] = reginput;
		return true;
	}
	else
		return false;
}

// regmatch - main matching routine
//
// Conceptually the strategy is simple:  check to see whether the current
// node matches, call self recursively to see whether the rest matches,
// and then act accordingly.  In practice we make some effort to avoid
// recursion, in particular by going through "ordinary" nodes (that don't
// need to know whether the rest of the match failed) by a loop instead of
// by recursion.

bool CRegExecutor::regmatch( char* prog )
{
	char* scan;	// Current node. 
	char* next;		// Next node. 

#ifdef _RE_DEBUG
	if (prog != NULL && regnarrate)
		_ftprintf(stderr, "%s(\n", regprop(prog));
#endif
	for (scan = prog; scan != NULL; scan = next)
	{
#ifdef _RE_DEBUG
		if (regnarrate)
			_ftprintf(stderr, "%s...\n", regprop(scan));
#endif
		next = regnext(scan);

		switch (OP(scan))
		{
			case BOL:
				if (reginput != regbol)
					return false;
				break;
			case EOL:
				if (*reginput != '\0')
					return false;
				break;
			case WORDA:
				/* Must be looking at a letter, digit, or _ */
				if ((!isalnum(*reginput)) && *reginput != '_')
					return(0);
				/* Prev must be BOL or nonword */
				if (reginput > regbol &&
					(isalnum(reginput[-1]) || reginput[-1] == '_'))
					return(0);
				break;
			case WORDZ:
				/* Must be looking at non letter, digit, or _ */
				if (isalnum(*reginput) || *reginput == '_')
					return(0);
				/* We don't care what the previous char was */
				break;
			case ANY:
				if (*reginput == '\0')
					return false;
				reginput++;
				break;
			case EXACTLY:
			{
				size_t len;
				char* const opnd = OPERAND(scan);

				// Inline the first character, for speed. 
				if (*opnd != *reginput)
					return false;
				len = strlen(opnd);
				if (len > 1 && strncmp(opnd, reginput, len) != 0)
					return false;
				reginput += len;

				break;
			}
			case ANYOF:
				if (*reginput == '\0' ||
						strchr(OPERAND(scan), *reginput) == NULL)
					return false;
				reginput++;
				break;
			case ANYBUT:
				if (*reginput == '\0' ||
						strchr(OPERAND(scan), *reginput) != NULL)
					return false;
				reginput++;
				break;
			case NOTHING:
				break;
			case BACK:
				break;
			case OPEN+1: case OPEN+2: case OPEN+3:
			case OPEN+4: case OPEN+5: case OPEN+6:
			case OPEN+7: case OPEN+8: case OPEN+9:
			{
				const int no = OP(scan) - OPEN;
				char* const input = reginput;

				if (regmatch(next))
				{
					// Don't set startp if some later
					// invocation of the same parentheses
					// already has.
					 
					if (regstartp[no] == NULL)
						regstartp[no] = input;
					return true;
				}
				else
					return false;
				break;
			}
			case CLOSE+1: case CLOSE+2: case CLOSE+3:
			case CLOSE+4: case CLOSE+5: case CLOSE+6:
			case CLOSE+7: case CLOSE+8: case CLOSE+9:
			{
				const int no = OP(scan) - CLOSE;
				char* const input = reginput;

				if (regmatch(next))
				{
					// Don't set endp if some later
					// invocation of the same parentheses
					// already has.
					 
					if (regendp[no] == NULL)
						regendp[no] = input;
					return true;
				}
				else
					return false;
				break;
			}
			case BRANCH:
			{
				char* const save = reginput;

				if (OP(next) != BRANCH)		// No choice. 
					next = OPERAND(scan);	// Avoid recursion. 
				else
				{
					while (OP(scan) == BRANCH)
					{
						if (regmatch(OPERAND(scan)))
							return true;
						reginput = save;
						scan = regnext(scan);
					}
					return false;
					// NOTREACHED 
				}
				break;
			}
			case STAR: case PLUS:
			{
				const char nextch = (OP(next) == EXACTLY) ? *OPERAND(next) : '\0';
				size_t no;
				char* const save = reginput;
				const size_t min = (OP(scan) == STAR) ? 0 : 1;

				for (no = regrepeat(OPERAND(scan)) + 1; no > min; no--)
				{
					reginput = save + no - 1;
					// If it could work, try it. 
					if (nextch == '\0' || *reginput == nextch)
						if (regmatch(next))
							return true;
				}
				return false;
				break;
			}
			case END:
				return true;	// Success! 
				break;
			default:
				regerror( REGERR_CORRUPTION );
				return false;
				break;
		}
	}

	// We get here only if there's trouble -- normally "case END" is
	// the terminating point.

	regerror( REGERR_CORRUPTED_POINTERS );
	return false;
}

// regrepeat - report how many times something simple would match

size_t CRegExecutor::regrepeat( char* node )
{
	size_t count;
	char* scan;
	char ch;

	switch (OP(node))
	{
		case ANY:
			return(strlen(reginput));
			break;
		case EXACTLY:
			ch = *OPERAND(node);
			count = 0;
			for (scan = reginput; *scan == ch; scan++)
				count++;
			return(count);
			break;
		case ANYOF:
			return(strspn(reginput, OPERAND(node)));
			break;
		case ANYBUT:
			return(strcspn(reginput, OPERAND(node)));
			break;
		default:		// Oh dear.  Called inappropriately. 
			regerror( REGERR_BAD_REGREPEAT );
			return(0);	// Best compromise. 
			break;
	}
	// NOTREACHED 
}

#ifdef _RE_DEBUG

// regdump - dump a regexp onto stdout in vaguely comprehensible form

void regexp_internal::regdump()
{
	char* s;
	char op = EXACTLY;	// Arbitrary non-END op. 
	char* next;

	s = strinc(program);
	while (op != END)
	{	// While that wasn't END last time... 
		op = OP(s);
		_tprintf("%2d%s", s-program, regprop(s));	// Where, what. 
		next = regnext(s);
		if (next == NULL)		// Next ptr. 
			_tprintf("(0)");
		else 
			_tprintf("(%d)", (s-program)+(next-s));
		s += 3;
		if (op == ANYOF || op == ANYBUT || op == EXACTLY)
		{
			// Literal string, where present. 
			while (*s != '\0')
			{
				_putchar(*s);
				s = strinc(s);
			}
			s = strinc(s);
		}
		_putchar('\n');
	}

	// Header fields of interest. 
	if (regstart != '\0')
		_tprintf("start `%c' ", regstart);
	if (reganch)
		_tprintf("anchored ");
	if (regmust != NULL)
		_tprintf("must have \"%s\"", regmust);
	_tprintf("\n");
}

// regprop - printable representation of opcode

#define OUTPUT(s) case s: p = s; break
char* CRegProgramAccessor::regprop( char* op )
{
	char* p;
	static char buf[50];

	(void) strcpy(buf, ":");

	switch (OP(op))
	{
	OUTPUT( BOL );
	OUTPUT( EOL );
	OUTPUT( ANY );
	OUTPUT( ANYOF );
	OUTPUT( ANYBUT );
	OUTPUT( BRANCH );
	OUTPUT( EXACTLY );
	OUTPUT( NOTHING );
	OUTPUT( BACK );
	OUTPUT( END );
	OUTPUT( STAR );
	OUTPUT( PLUS );
	OUTPUT( WORDA );
	OUTPUT( WORDZ );
	case OPEN+1: case OPEN+2: case OPEN+3:
	case OPEN+4: case OPEN+5: case OPEN+6:
	case OPEN+7: case OPEN+8: case OPEN+9:
		_stprintf(buf+strlen(buf), "OPEN%d", OP(op)-OPEN);
		p = NULL;
		break;
	case CLOSE+1: case CLOSE+2: case CLOSE+3:
	case CLOSE+4: case CLOSE+5: case CLOSE+6:
	case CLOSE+7: case CLOSE+8: case CLOSE+9:
		_stprintf(buf+strlen(buf), "CLOSE%d", OP(op)-CLOSE);
		p = NULL;
		break;
	default:
		regerror( REGERR_CORRUPTED_OPCODE );
		break;
	}
	if (p != NULL)
		(void) strcat(buf, p);
	return(buf);
}
#endif

///////////////////////////////////////////////////////////////////////////////

// constructor
Regexp::Regexp()
	: str(0),
	  rc(0)
{ }

Regexp::Regexp( const char* exp, bool iCase )
	: str( 0 ),
	  rc( new regexp_internal( exp, iCase ) )
{ }

Regexp::Regexp( const std::string& exp, bool iCase )
	: str( 0 ),
	  rc( new regexp_internal( exp.c_str(), iCase ) )
{ }

Regexp::Regexp( const Regexp &r )
  : str(r.str),
    m_szError(r.m_szError),
    rc( r.rc )
{
	if ( rc )
		rc->count++;
}

const Regexp & Regexp::operator=( const Regexp & r )
{
	if ( this != &r )
	{
		if ( rc && rc->count-- == 0 )
			delete rc;

		rc = r.rc;
		if ( rc )
			rc->count++;

		str = r.str;
		m_szError = r.m_szError;
	}	
	return *this;
}

Regexp::~Regexp()
{
	if ( rc && rc->count-- == 0 )
		delete rc;
}

bool Regexp::Match( const char * s )
{
	ClearErrorString();
	str = s;
	bool ret = false;
	if ( rc )
	{
		// copy on write !

		if ( rc->count )
		{
			rc->count--;
			rc = rc->getCopy();
		}

		ret = rc->regexec( (char*) s );
		int i = 0;
		if ( ret )
			for ( i = 0; i < Regexp::NSUBEXP && rc->startp[i] ; i++ )
				;
		rc->numSubs = i - 1;
	}
	else
		m_szError = CRegErrorHandler::FindErr( REGERR_NO_REGEXP );
	return ret;
}

std::string Regexp::GetReplaceString( char* source ) const
{
	ClearErrorString();
	if ( rc )
		return rc->GetReplaceString( source );
	else
		m_szError = CRegErrorHandler::FindErr( REGERR_NO_REGEXP );
	return "";
}

int Regexp::SubStrings() const
{
	ClearErrorString();
	int ret = -1;
	if ( rc )
		ret = rc->numSubs;
	else
		m_szError = CRegErrorHandler::FindErr( REGERR_NO_REGEXP );
	return ret;
}

int Regexp::SubStart( unsigned int i ) const
{
	ClearErrorString();
	int ret = -1;
	if ( rc )
		ret = rc->startp[safeIndex(i)] - str;
	else
		m_szError = CRegErrorHandler::FindErr( REGERR_NO_REGEXP );
	return ret;
}

int Regexp::SubLength( unsigned int i ) const
{
	ClearErrorString();
	int ret = -1;
	if ( rc )
	{
		i = safeIndex(i);
		ret = rc->endp[i] - rc->startp[i];
	}
	else
		m_szError = CRegErrorHandler::FindErr( REGERR_NO_REGEXP );
	return ret;
}

bool Regexp::CompiledOK() const
{
	return rc ? rc->Status() : false;
}

#ifdef _RE_DEBUG
void Regexp::Dump()
{
	if ( rc )
		rc->regdump();
#if defined( _DEBUG )
	else
		TRACE0( "No regexp to dump out\n" );
#endif
}
#endif

int Regexp::safeIndex( unsigned int i ) const
{
	return i < Regexp::NSUBEXP ? i : Regexp::NSUBEXP;
}

const std::string Regexp::operator[]( unsigned int i ) const
{
	ClearErrorString();
	assert( rc );
	if ( rc )
	{
		int len = SubLength(i);
		std::string buffer (rc->startp[i], len);
		return buffer;
	}
	else
	{
		m_szError = CRegErrorHandler::FindErr( REGERR_NO_REGEXP );
		return "";
	}
}

void regexp_internal::ignoreCase( const char * in, char * out )
{
	// copy in to out making every top level character a [Aa] set
	bool inRange = 0;
	while( *in )
	{
		if ( *in == '[' )
			inRange = 1;
		if ( *in == ']' )
			inRange = 0;
		if ( ! inRange && isalpha( *in ) )
		{
			*out++ = '[';
			*out++ = (char)toupper( *in );
			*out++ = (char)tolower( *in );
			*out++ = ']';
		}
		else
			*out++ = *in;
		in++;
	}
	*out = 0;
}

// GetReplaceString	- Converts a replace expression to a string
//					- perform substitutions after a regexp match
// Returns			- The resultant string
std::string regexp_internal::GetReplaceString( const char* sReplaceExp ) const
{
	std::string szEmpty;

	char *src = (char *)sReplaceExp;
	char *buf = 0;
	char c;
	int no;
	size_t len;

	if( sReplaceExp == NULL )
	{
		regerror( REGERR_NULL_TO_REGSUB  );
		return szEmpty;
	}
	if ( *program != MAGIC)
	{
		regerror( REGERR_DAMAGED_REGEXP_REGSUB );
		return szEmpty;
	}

	// First compute the length of the string
	int replacelen = 0;
	while ((c = *src++) != '\0') 
	{
		if (c == '&')
			no = 0;
		else if (c == '\\' && isdigit(*src))
			no = *src++ - '0';
		else
			no = -1;

		if (no < 0) 
		{	
			// Ordinary character. 
			if (c == '\\' && (*src == '\\' || *src == '&'))
				c = *src++;
			replacelen++;
		} 
		else if (startp[no] != NULL && endp[no] != NULL &&
					endp[no] > startp[no]) 
		{
			// Get tagged expression
			len = endp[no] - startp[no];
			replacelen += len;
		}
	}

	std::string szReplace (replacelen, '\0');

	// Now we can create the string
	src = (char *)sReplaceExp;
	while ((c = *src++) != '\0') 
	{
		if (c == '&')
			no = 0;
		else if (c == '\\' && isdigit(*src))
			no = *src++ - '0';
		else
			no = -1;

		if (no < 0) 
		{	
			// Ordinary character. 
			if (c == '\\' && (*src == '\\' || *src == '&'))
				c = *src++;
			szReplace.push_back(c);
		} 
		else if (startp[no] != NULL && endp[no] != NULL &&
					endp[no] > startp[no]) 
		{
			// Get tagged expression
			len = endp[no] - startp[no];
			strncpy(buf, startp[no], len);
			buf += len;
			if (len != 0 && *(buf-1) == '\0')
			{	/* strncpy hit NUL. */
				regerror( REGERR_DAMAGED_MATCH_STRING );
				return szEmpty;
			}
		}
	}

	return szReplace;
}

std::string Regexp::GetErrorString() const
{
	// make sure that if status == 0 that we have an error string
	assert( ( ! CompiledOK() ) ? ( rc ? rc->GetErrorString() : m_szError).size() != 0 : 1 );
	return rc ? rc->GetErrorString() : m_szError ;
}

void Regexp::ClearErrorString() const
{
	if ( rc )
		rc->ClearErrorString();
	m_szError.clear();
}

