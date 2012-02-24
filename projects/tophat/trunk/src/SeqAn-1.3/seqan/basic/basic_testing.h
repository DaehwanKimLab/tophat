// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// The SeqAn testing infrastructure.  Based on ideas from the OpenMS
// "ClassTest.h".
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef SEQAN_BASIC_BASIC_TESTING_H_
#define SEQAN_BASIC_BASIC_TESTING_H_

#include <iostream>  // stdout, stderr
#include <cstring>   // strrpos
#include <cstdlib>   // exit()
#include <cstdarg>   // va_start, va_list, va_end
#include <set>
#include <vector>
#include <string>

#ifdef PLATFORM_WINDOWS
#include <Windows.h>  // DeleteFile()
#else  // #ifdef PLATFORM_WINDOWS
#include <unistd.h>  // unlink()
#endif  // #ifdef PLATFORM_WINDOWS

// SeqAn's has three global debug/testing levels: testing, debug and
// release.  Depending on the level, the SEQAN_ASSERT_* and
// SEQAN_CHECKPOINT macros will be enabled.
//
// Note that this is independent of the <cassert> assertions and
// NDEBUG being defined.
//
// The levels are enabled by the values of the macros
// SEQAN_ENABLE_TESTING and SEQAN_ENABLE_DEBUG.  By setting a macro to
// 0, one disables the level and by setting the macro to 1, one
// enables a level.  Enabling testing also enables debug, overriding a
// value of 0 for SEQAN_ENABLE_DEBUG.
//
// If the level is release (both the macros for debug and testing are
// 0), the assertions will be disabled.  If the level is debug then
// the assertions will be enabled.  If the level is testing then the
// checkpoint macros will also be enabled.
//
// The default is to enable debugging but disable testing.
//
// You can print the current level using the function seqan::printDebugLevel().

// Set default for SEQAN_ENABLE_TESTING.
#ifndef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 0
#endif  // #ifndef SEQAN_ENABLE_TESTING

// Set default for SEQAN_ENABLE_DEBUG.
#ifndef SEQAN_ENABLE_DEBUG
#define SEQAN_ENABLE_DEBUG 1
#endif  // #ifndef SEQAN_ENABLE_DEBUG

// Force-enable debugging if testing is enabled.
#if SEQAN_ENABLE_TESTING
#undef SEQAN_ENABLE_DEBUG
#define SEQAN_ENABLE_DEBUG 1
#endif  // #if SEQAN_ENABLE_TESTING

// Allow disabling checkpoints independent of testing.
#ifndef SEQAN_ENABLE_CHECKPOINTS
#define SEQAN_ENABLE_CHECKPOINTS SEQAN_ENABLE_TESTING
#endif  // #ifndef SEQAN_ENABLE_CHECKPOINTS

namespace seqan {

// SEQAN_CXX_FLAGS_ contains the compiler flags, SEQAN_CXX_FLAGS is a string
// literal with this value.
#if !defined(SEQAN_CXX_FLAGS_)
#define SEQAN_CXX_FLAGS_ SEQAN_CXX_FLAGS_NOT_SET
#endif //  !defined(SEQAN_CXX_FLAGS__)
#define SEQAN_MKSTRING_(str) # str
#define SEQAN_MKSTRING(str) SEQAN_MKSTRING_(str)
#define SEQAN_CXX_FLAGS SEQAN_MKSTRING(SEQAN_CXX_FLAGS_)
//#undef SEQAN_MKSTRING
//#undef SEQAN_MKSTRING_


/**
.Function.printDebugLevel:
..cat:Miscellaneous:
..summary:Print the current SeqAn debug level and the compiler flags to the given stream.
..signature:printDebugLevel(stream)
..param.stream:The stream to print to, e.g. $std::cout$.
..include:seqan/basic.h
 */
template <typename TStream>
void printDebugLevel(TStream &stream) {
    stream << "SEQAN_ENABLE_DEBUG == " << SEQAN_ENABLE_DEBUG << std::endl;
    stream << "SEQAN_ENABLE_TESTING == " << SEQAN_ENABLE_TESTING << std::endl;
    stream << "SEQAN_ENABLE_CHECKPOINTS == " << SEQAN_ENABLE_CHECKPOINTS << std::endl;
    stream << "SEQAN_CXX_FLAGS == \"" << SEQAN_CXX_FLAGS << "\"" << std::endl;
}


// Namespace for the testing infrastructure.
//
// This namespace contains the variables and functions that are used
// in the macros below to perform the tests.
namespace ClassTest {
    // Raised when an assertion fails in test mode.
    struct AssertionFailedException {};

    // Container for static global data for the tests.
    struct StaticData {
        // Number of tests that were run.
        static int &testCount() {
            static int result = 0;
            return result;
        }

        // Number of errors that occured.
        static int &errorCount() {
            static int result = 0;
            return result;
        }

        // Number of skipped tests.
        static int &skippedCount() {
            static int result = 0;
            return result;
        }

        // Flag whether there was an error in this test.
        static bool &thisTestOk() {
            static bool result = 0;
            return result;
        }

        // Flag whether this test was skipped.
        static bool &thisTestSkipped() {
            static bool result = 0;
            return result;
        }

        // Name of the current test.
        static const char *&currentTestName() {
            const char *defaultValue = "";
            static const char *result = const_cast<char*>(defaultValue);
            return result;
        }

        // Base path to the binary.  Extrapolated from __FILE__.
        static char *&basePath() {
            const char *defaultValue = ".";
            static char *result = const_cast<char*>(defaultValue);
            return result;
        }

        // Base path to the "projects" directory, extrapolated from
        // __FILE__.
        static char *&pathToProjects() {
            const char *defaultValue = ".";
            static char *result = const_cast<char*>(defaultValue);
            return result;
        }

        // Total number of checkpoints in header file.
        static int &totalCheckPointCount() {
            static int result = 0;
            return result;
        }

        // Total number of checkpoints found in binary files.
        static int &foundCheckPointCount() {
            static int result = 0;
            return result;
        }

        // Names of temporary files as returned by tempFileName.  This
        // global state is used to remove any existing such files
        // after completing the testsuite.
        static ::std::vector<std::string> & tempFileNames() {
            static ::std::vector<std::string> filenames;
            return filenames;
        }
    };

// Open a temporary file, unlink it, return posix handle.  Note: This has not been tested yet.
// TODO(holtgrew): Not used yet and Windows code does not work.
/*
inline
int openTempFile() {
#ifdef PLATFORM_WINDOWS
    char * fileName = _tempnam(NULL, "SQN");
    if (!fileName) {
        ::std::cerr << "Cannot create a unique temporary filename" << ::std::endl;
        exit(1);
    }
    int result = open(fileName, _O_RDWR | OPEN_TEMPORARY);
    free(fileName);
    return result;
#else  // A Unix...
    char filenameBuffer[100];
    strcpy(filenameBuffer, "/tmp/SEQANXXXXXXXXXX");
    int result = mkstemp(filenameBuffer);
    unlink(filenameBuffer);
    return result;
#endif  // ifdef PLATFORM_WINDOWS
}
*/

// Return the path to a temporary file, in a static buffer in this
// function.  This is not thread safe!
inline
const char *tempFileName() {
    static char fileNameBuffer[100];
#ifdef PLATFORM_WINDOWS_VS
    char * fileName = tempnam(NULL, "SEQAN.");
    if (!fileName) {
        ::std::cerr << "Cannot create a unique temporary filename" << ::std::endl;
        exit(1);
    }
    strcpy(fileNameBuffer, fileName);
    free(fileName);
    StaticData::tempFileNames().push_back(fileNameBuffer);
    return fileNameBuffer;
#else  // ifdef PLATFORM_WINDOWS_VS
    strcpy(fileNameBuffer, "/tmp/SEQAN.XXXXXXXXXXXXXXXXXXXX");
#ifdef PLATFORM_WINDOWS_MINGW
    // There is no mkstemp in MinGW but it does not complain about tmpnam.
    tmpnam(fileNameBuffer);
#else  // ifdef PLATFORM_WINDOWS_MINGW
    mkstemp(fileNameBuffer);
    unlink(fileNameBuffer);
#endif  // #ifdef PLATFORM_WINDOWS_MINGW
    StaticData::tempFileNames().push_back(fileNameBuffer);
    return fileNameBuffer;
#endif  // ifdef PLATFORM_WINDOWS_VS
}

    // Initialize the testing infrastructure.
    //
    // Used through SEQAN_BEGIN_TESTSUITE(test_name)
    inline
    void beginTestSuite(const char *testSuiteName, const char *argv0) {
        // First things first: Print the current debug level.
        printDebugLevel(std::cout);
        (void)testSuiteName;
        StaticData::testCount() = 0;
        StaticData::skippedCount() = 0;
        StaticData::errorCount() = 0;
        StaticData::totalCheckPointCount() = 0;
        StaticData::foundCheckPointCount() = 0;
        // Get path to argv0.
        const char *end = argv0;
#ifdef PLATFORM_WINDOWS
        const char pathSeparator = '\\';
#else  // PLATFORM_WINDOWS
        const char pathSeparator = '/';
#endif  // PLATFORM_WINDOWS
        for (const char *ptr = strchr(argv0, pathSeparator); ptr != 0; ptr = strchr(ptr+1, pathSeparator))
            end = ptr;
        int rpos = end - argv0;
        if (rpos <= 0) {
            StaticData::basePath() = new char[1];
            strcpy(StaticData::basePath(), ".");
        } else {
            int len = rpos;
            StaticData::basePath() = new char[len];
            strncpy(StaticData::basePath(), argv0, len);
        }
        // Get path to projects.
        const char *file = __FILE__;
        int pos = -1;
        for (size_t i = 0; i < strlen(file) - strlen("projects"); ++i) {
            if (strncmp(file + i, "projects", strlen("projects")) == 0) {
                pos = i;
            }
        }
        if (pos == -1) {
            std::cerr << "Could not extrapolate path to projects from __FILE__ == \""
                      << __FILE__ << "\"" << std::endl;
            exit(1);
        }
        StaticData::pathToProjects() = new char[pos];
        strncpy(StaticData::pathToProjects(), file, pos);
        StaticData::pathToProjects()[pos-1] = '\0';
    }

    // Run test suite finalization.
    //
    // Used through SEQAN_END_TESTSUITE
    //
    // Prints a bottom banner with the error count and returns the
    // program's return code.
    inline
    int endTestSuite() {
        delete[] StaticData::basePath();
        delete[] StaticData::pathToProjects();

        std::cout << "**************************************" << std::endl;
        std::cout << " Total Check Points : " << StaticData::totalCheckPointCount() << std::endl;
        std::cout << " Found Check Points : " << StaticData::foundCheckPointCount() << std::endl;
        std::cout << " Lost Check Points  : " << StaticData::totalCheckPointCount() - StaticData::foundCheckPointCount() << std::endl;
        std::cout << "--------------------------------------" << std::endl;
        std::cout << " Total Tests: " << StaticData::testCount() << std::endl;
        std::cout << " Skipped:     " << StaticData::skippedCount() << std::endl;
        std::cout << " Errors:      " << StaticData::errorCount() << std::endl;
        std::cout << "**************************************" << std::endl;
        if (StaticData::errorCount() != 0)
            return 1;
        // TODO(holtgrew): Re-enable that all check points have to be found for the test to return 1;
        /*
        if (StaticData::totalCheckPointCount() != StaticData::foundCheckPointCount())
            return 1;
        */
        // Delete all temporary files that still exist.
        for (unsigned i = 0; i < StaticData::tempFileNames().size(); ++i) {
#ifdef PLATFORM_WINDOWS
            DeleteFile(StaticData::tempFileNames()[i].c_str());
#else  // #ifdef PLATFORM_WINDOWS
            unlink(StaticData::tempFileNames()[i].c_str());
#endif  // #ifdef PLATFORM_WINDOWS
        }

        return 0;
    }

    // Run test initialization.
    inline
    void beginTest(const char *testName) {
        StaticData::currentTestName() = testName;
        StaticData::thisTestOk() = true;
        StaticData::thisTestSkipped() = false;
        StaticData::testCount() += 1;
    }

    // Run test finalization.
    inline
    void endTest() {
        if (StaticData::thisTestSkipped()) {
            std::cout << StaticData::currentTestName() << " SKIPPED" << std::endl;
        } else if (StaticData::thisTestOk()) {
            std::cout << StaticData::currentTestName() << " OK" << std::endl;
        } else {
            std::cerr << StaticData::currentTestName() << " FAILED" << std::endl;
        }
    }

    // Marks the current test as "skipped".
    inline
    void skipCurrentTest() {
        StaticData::thisTestSkipped() = true;
        StaticData::skippedCount() += 1;
    }

    // Called by the macro SEQAN_ASSERT_FAIL.
    inline void forceFail(const char *file, int line,
                          const char *comment, ...) {
        StaticData::errorCount() += 1;
        std::cerr << file << ":" << line << " FAILED! ";
        if (comment) {
            std::cerr << " (";
            va_list args;
            va_start(args, comment);
            vfprintf(stderr, comment, args);
            va_end(args);
            std::cerr << ")";
        }
        std::cerr << std::endl;
    }

    // Similar to forceFail above, but accepting a va_list parameter.
    inline void vforceFail(const char *file, int line,
                           const char *comment, va_list argp) {
        StaticData::errorCount() += 1;
        std::cerr << file << ":" << line << " FAILED! ";
        if (comment) {
            std::cerr << " (";
            vfprintf(stderr, comment, argp);
            std::cerr << ")";
        }
        std::cerr << std::endl;
    }

    // Same as forceFail above, but with comment set to 0.
    inline void forceFail(const char *file, int line) {
        forceFail(file, line, 0);
    }

    // Called by the macro SEQAN_ASSERT_EQ.
    //
    // Tests that the given two value are equal.  Returns true iff the
    // two values are equal.
    template <typename T1, typename T2>
    bool testEqual(const char *file, int line,
                   const T1 &value1, const char *expression1,
                   const T2 &value2, const char *expression2,
                   const char *comment, ...) {
        if (!(value1 == value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " == " << expression2 << " was: " << value1
                      << " != " << value2;
            if (comment) {
                std::cerr << " (";
                va_list args;
                va_start(args, comment);
                vfprintf(stderr, comment, args);
                va_end(args);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Similar to testEqual above, but accepts a va_list instead of variadic
    // parameters.
    template <typename T1, typename T2>
    bool vtestEqual(const char *file, int line,
                    const T1 &value1, const char *expression1,
                    const T2 &value2, const char *expression2,
                    const char *comment, va_list argp) {
        if (!(value1 == value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " == " << expression2 << " was: " << value1
                      << " != " << value2;
            if (comment) {
                std::cerr << " (";
                vfprintf(stderr, comment, argp);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Same as testEqual above, but with comment set to 0.
    template <typename T1, typename T2>
    bool testEqual(const char *file, int line,
                   const T1 &value1, const char *expression1,
                   const T2 &value2, const char *expression2) {
        return testEqual(file, line, value1, expression1, value2, expression2, 0);
    }



    // Called by the macro SEQAN_ASSERT_IN_DELTA.
    //
    // Tests that the given two value are equal.  Returns true iff the
    // two values are equal.
    template <typename T1, typename T2, typename T3>
    bool testInDelta(const char *file, int line,
                     const T1 &value1, const char *expression1,
                     const T2 &value2, const char *expression2,
                     const T3 &value3, const char *expression3,
                     const char *comment, ...) {
        if (!(value1 >= value2 - value3 && value1 <= value2 + value3)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " in [" << expression2 << " - " << expression3
                      << ", " << expression2 << " + " << expression3 << "] was: " << value1
                      << " not in [" << value2 - value3 << ", " << value2 + value3 << "]";
            if (comment) {
                std::cerr << " (";
                va_list args;
                va_start(args, comment);
                vfprintf(stderr, comment, args);
                va_end(args);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Similar to testInDelta above, but accepts a va_list instead of variadic
    // parameters.
    template <typename T1, typename T2, typename T3>
    bool vtestInDelta(const char *file, int line,
                      const T1 &value1, const char *expression1,
                      const T2 &value2, const char *expression2,
                      const T3 &value3, const char *expression3,
                      const char *comment, va_list argp) {
        if (!(value1 >= value2 - value3 && value1 <= value2 + value3)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " in [" << expression2 << " - " << expression3
                      << ", " << expression2 << " + " << expression3 << "] was: " << value1
                      << " not in [" << value2 - value3 << ", " << value2 + value3 << "]";
            if (comment) {
                std::cerr << " (";
                vfprintf(stderr, comment, argp);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Same as testInDelta above, but with comment set to 0.
    template <typename T1, typename T2, typename T3>
    bool testInDelta(const char *file, int line,
                     const T1 &value1, const char *expression1,
                     const T2 &value2, const char *expression2,
                     const T3 &value3, const char *expression3) {
        return testInDelta(file, line, value1, expression1, value2, expression2, value3, expression3, 0);
    }


    // Called by the macro SEQAN_ASSERT_NEQ.
    //
    // Tests that the given two value are not equal.  Returns true iff
    // the two values are equal.
    template <typename T1, typename T2>
    bool testNotEqual(const char *file, int line,
                      const T1 &value1, const char *expression1,
                      const T2 &value2, const char *expression2,
                      const char *comment, ...) {
        if (!(value1 != value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " != " << expression2 << " was: " << value1
                      << " == " << value2;
            if (comment) {
                std::cerr << " (";
                va_list args;
                va_start(args, comment);
                vfprintf(stderr, comment, args);
                va_end(args);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }

    // Similar to testNotEqual above, but accepts a va_list instead of variadic
    // parameters.
    template <typename T1, typename T2>
    bool vtestNotEqual(const char *file, int line,
                       const T1 &value1, const char *expression1,
                       const T2 &value2, const char *expression2,
                       const char *comment, va_list argp) {
        if (!(value1 != value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " != " << expression2 << " was: " << value1
                      << " == " << value2;
            if (comment) {
                std::cerr << " (";
                vfprintf(stderr, comment, argp);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Same as testNotEqual above, but with comment set to 0.
    template <typename T1, typename T2>
    bool testNotEqual(const char *file, int line,
                      const T1 &value1, const char *expression1,
                      const T2 &value2, const char *expression2) {
        return testNotEqual(file, line, value1, expression1, value2, expression2, 0);
    }


    // Called by the macro SEQAN_ASSERT_GEQ.
    //
    // Tests that the first value is greater than or equal to the
    // second one.  Returns true iff the test yields true.
    template <typename T1, typename T2>
    bool testGeq(const char *file, int line,
                 const T1 &value1, const char *expression1,
                 const T2 &value2, const char *expression2,
                 const char *comment, ...) {
        if (!(value1 >= value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " >= " << expression2 << " was: " << value1
                      << " < " << value2;
            if (comment) {
                std::cerr << " (";
                va_list args;
                va_start(args, comment);
                vfprintf(stderr, comment, args);
                va_end(args);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Similar to testGeq above, but accepts a va_list instead of variadic
    // parameters.
    template <typename T1, typename T2>
    bool vtestGeq(const char *file, int line,
                  const T1 &value1, const char *expression1,
                  const T2 &value2, const char *expression2,
                  const char *comment, va_list argp) {
        if (!(value1 >= value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " >= " << expression2 << " was: " << value1
                      << " < " << value2;
            if (comment) {
                std::cerr << " (";
                vfprintf(stderr, comment, argp);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }

    // Same as testGeq above, but with comment set to 0.
    template <typename T1, typename T2>
    bool testGeq(const char *file, int line,
                 const T1 &value1, const char *expression1,
                 const T2 &value2, const char *expression2) {
        return testGeq(file, line, value1, expression1, value2, expression2, 0);
    }


    // Called by the macro SEQAN_ASSERT_GT.
    //
    // Tests that the first value is greater than the second one.
    // Returns true iff the test yields true.
    template <typename T1, typename T2>
    bool testGt(const char *file, int line,
                const T1 &value1, const char *expression1,
                const T2 &value2, const char *expression2,
                const char *comment, ...) {
        if (!(value1 > value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " > " << expression2 << " was: " << value1
                      << " <= " << value2;
            if (comment) {
                std::cerr << " (";
                va_list args;
                va_start(args, comment);
                vfprintf(stderr, comment, args);
                va_end(args);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Similar to testGt above, but accepts a va_list instead of variadic
    // parameters.
    template <typename T1, typename T2>
    bool vtestGt(const char *file, int line,
                 const T1 &value1, const char *expression1,
                 const T2 &value2, const char *expression2,
                 const char *comment, va_list argp) {
        if (!(value1 > value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " > " << expression2 << " was: " << value1
                      << " <= " << value2;
            if (comment) {
                std::cerr << " (";
                vfprintf(stderr, comment, argp);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }

                
    // Same as testGt above, but with comment set to 0.
    template <typename T1, typename T2>
    bool testGt(const char *file, int line,
                const T1 &value1, const char *expression1,
                const T2 &value2, const char *expression2) {
        return testGt(file, line, value1, expression1, value2, expression2, 0);
    }


    // Called by the macro SEQAN_ASSERT_LEQ.
    //
    // Tests that the first value is less than or equal to the second
    // one.  Returns true iff the test yields true.
    template <typename T1, typename T2>
    bool testLeq(const char *file, int line,
                 const T1 &value1, const char *expression1,
                 const T2 &value2, const char *expression2,
                 const char *comment, ...) {
        if (!(value1 <= value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " <= " << expression2 << " was: " << value1
                      << " > " << value2;
            if (comment) {
                std::cerr << " (";
                va_list args;
                va_start(args, comment);
                vfprintf(stderr, comment, args);
                va_end(args);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Similar to testLeq above, but accepts a va_list instead of variadic
    // parameters.
    template <typename T1, typename T2>
    bool vtestLeq(const char *file, int line,
                  const T1 &value1, const char *expression1,
                  const T2 &value2, const char *expression2,
                  const char *comment, va_list argp) {
        if (!(value1 <= value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " <= " << expression2 << " was: " << value1
                      << " > " << value2;
            if (comment) {
                std::cerr << " (";
                vfprintf(stderr, comment, argp);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Same as testLeq above, but with comment set to 0.
    template <typename T1, typename T2>
    bool testLeq(const char *file, int line,
                 const T1 &value1, const char *expression1,
                 const T2 &value2, const char *expression2) {
        return testLeq(file, line, value1, expression1, value2, expression2, 0);
    }


    // Called by the macro SEQAN_ASSERT_LT.
    //
    // Tests that the first value is greater than the second one.
    // Returns true iff the test yields true.
    template <typename T1, typename T2>
    bool testLt(const char *file, int line,
                const T1 &value1, const char *expression1,
                const T2 &value2, const char *expression2,
                const char *comment, ...) {
        if (!(value1 < value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " < " << expression2 << " was: " << value1
                      << " >= " << value2;
            if (comment) {
                std::cerr << " (";
                va_list args;
                va_start(args, comment);
                vfprintf(stderr, comment, args);
                va_end(args);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Similar to testLt above, but accepts a va_list instead of variadic
    // parameters.
    template <typename T1, typename T2>
    bool vtestLt(const char *file, int line,
                 const T1 &value1, const char *expression1,
                 const T2 &value2, const char *expression2,
                 const char *comment, va_list argp) {
        if (!(value1 < value2)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " < " << expression2 << " was: " << value1
                      << " >= " << value2;
            if (comment) {
                std::cerr << " (";
                vfprintf(stderr, comment, argp);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Same as testLt above, but comment is 0.
    template <typename T1, typename T2>
    bool testLt(const char *file, int line,
                const T1 &value1, const char *expression1,
                const T2 &value2, const char *expression2) {
        return testLt(file, line, value1, expression1, value2, expression2, 0);
    }


    // Called by the macro SEQAN_ASSERT.
    //
    // Test that the given argument evaluates to true.
    template <typename T>
    bool testTrue(const char *file, int line,
                  const T &value_, const char *expression_,
                  const char *comment, ...) {
        if (!(value_)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression_ << " should be true but was " << (value_);
            if (comment) {
                std::cerr << " (";
                va_list args;
                va_start(args, comment);
                vfprintf(stderr, comment, args);
                va_end(args);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Similar to testTrue above, but accepts a va_list instead of variadic
    // parameters.
    template <typename T>
    bool vtestTrue(const char *file, int line,
                   const T &value_, const char *expression_,
                   const char *comment, va_list argp) {
        if (!(value_)) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression_ << " should be true but was " << (value_);
            if (comment) {
                std::cerr << " (";
                vfprintf(stderr, comment, argp);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Same as testTrue above, but comment will automatically be set to 0.
    template <typename T>
    bool testTrue(const char *file, int line,
                  const T &value_, const char *expression_)
    {
        return testTrue(file, line, value_, expression_, 0);
    }


    // Called by the macro SEQAN_ASSERT.
    //
    // Test that the given argument evaluates to false.
    template <typename T>
    bool testFalse(const char *file, int line,
                   const T &value_, const char *expression_,
                   const char *comment, ...) {
        if (value_) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression_ << " should be false but was " << (value_);
            if (comment) {
                std::cerr << " (";
                va_list args;
                va_start(args, comment);
                vfprintf(stderr, comment, args);
                va_end(args);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Similar to testFalse above, but accepts a va_list instead of variadic
    // parameters.
    template <typename T>
    bool vtestFalse(const char *file, int line,
                    const T &value_, const char *expression_,
                    const char *comment, va_list argp) {
        if (value_) {
            // Increase global error count.
            StaticData::thisTestOk() = false;
            StaticData::errorCount() += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression_ << " should be false but was " << (value_);
            if (comment) {
                std::cerr << " (";
                vfprintf(stderr, comment, argp);
                std::cerr << ")";
            }
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Same as testFalse above, but comment will automatically be set to 0.
    template <typename T>
    bool testFalse(const char *file, int line,
                   const T &value_, const char *expression_) {
        return testFalse(file, line, value_, expression_, 0);
    }

    // Represents a check point in a file.
    struct CheckPoint {
        // Path to the file.
        const char *file;
        // Line in the file.
        unsigned int line;

        // Less-than comparator for check points.
        bool operator<(const CheckPoint &other) const {
            int c = strcmp(file, other.file);
            if (c < 0)
                return true;
            if (c == 0 && line < other.line)
                return true;
            return false;
        }
    };

    // Wrapper for a set of check points.
    // TODO(holtgrew): Simply store the set?
    struct CheckPointStore {
        static ::std::set<CheckPoint> &data() {
            static ::std::set<CheckPoint> result;
            return result;
        }
    };

    // Puts the given check point into the CheckPointStore's data.
    inline bool
    registerCheckPoint(unsigned int line, const char *file) {
        const char *file_name = strrchr(file, '/');
        const char *file_name_2 = strrchr(file, '\\');
        if (file_name_2 > file_name)
            file_name = file_name_2;
        if (!file_name)
            file_name = file;
        else ++file_name;

        CheckPoint cp = {file_name, line};
        #ifdef _OMP
        #pragma omp critical
        #endif  // #ifdef _OMP
        CheckPointStore::data().insert(cp);
        return true;
    }

    // Test whether the given check point exists in the check point
    // store.
    inline void
    testCheckPoint(const char *file, unsigned int line) {
        StaticData::totalCheckPointCount() += 1;
        CheckPoint cp = {file, line};
        if (CheckPointStore::data().find(cp) == CheckPointStore::data().end()) {
            std::cerr << file << ":" << line << "  -- Check point lost."
                      << std::endl;
            return;
        }
        StaticData::foundCheckPointCount() += 1;
    }

    // Verify the check points for the given file.
    inline void
    verifyCheckPoints(const char *file) {
        char const* file_name = strrchr(file, '/');
        char const* file_name_2 = strrchr(file, '\\');
        if (file_name_2 > file_name) file_name = file_name_2;
        if (!file_name) file_name = file;
        else ++file_name;



        int len = strlen(StaticData::pathToProjects()) +
            strlen("/") + strlen(file) + 1;
        char *absolutePath = new char[len];
        absolutePath[0] = '\0';
        strcat(absolutePath, StaticData::pathToProjects());
        strcat(absolutePath, "/");
        strcat(absolutePath, file);

        FILE * fl = ::std::fopen(absolutePath, "r");
        delete[] absolutePath;
        if (!fl) {
            std::cerr << file << " -- verifyCheckPoints could not find this file." << std::endl;
        }
        unsigned int line_number = 1;
        char buf[1<<16];

        while (::std::fgets(buf, sizeof(buf), fl)) {
            if (::std::strstr(buf, "SEQAN_CHECKPOINT")) {
                testCheckPoint(file_name, line_number);
            }
            ++line_number;
        }

        ::std::fclose(fl);
    }

#if SEQAN_ENABLE_TESTING
    // If in testing mode then raise an AssertionFailedException.
    inline void fail() {
        StaticData::thisTestOk() = false;
        throw AssertionFailedException();
    }
#else
    // If not in testing mode then quit with an abort.
    inline void fail() {
        abort();
    }
#endif  // #if SEQAN_ENABLE_TESTING

}  // namespace ClassTest

// This macro expands to function header for one test.
#define SEQAN_DEFINE_TEST(test_name)                    \
    template <bool speed_up_dummy_to_prevent_compilation_of_unused_tests_> void SEQAN_TEST_ ## test_name ()

#if SEQAN_ENABLE_TESTING
// This macro expands to startup code for a test file.
#define SEQAN_BEGIN_TESTSUITE(suite_name)                       \
    int main(int argc, char **argv) {                           \
    (void) argc;                                                \
    ::seqan::ClassTest::beginTestSuite(#suite_name, argv[0]);


// This macro expands to shutdown code for a test file.
#define SEQAN_END_TESTSUITE                     \
    return ::seqan::ClassTest::endTestSuite();  \
}


// This macro expands to code to call a given test.
#define SEQAN_CALL_TEST(test_name)                                      \
    do {                                                                \
        ::seqan::ClassTest::beginTest(#test_name);                      \
        try {                                                           \
            SEQAN_TEST_ ## test_name<true>();                           \
        } catch(::seqan::ClassTest::AssertionFailedException e) {       \
            /* Swallow exception, go on with next test. */              \
            (void) e;  /* Get rid of unused variable warning. */        \
        }                                                               \
        ::seqan::ClassTest::endTest();                                  \
    } while (false)


// This macro returns from the current function and logs a "skipped"
// event for the current test.
#define SEQAN_SKIP_TEST                         \
    do {                                        \
        ::seqan::ClassTest::skipCurrentTest();  \
        return;                                 \
    } while (false)
#endif  // #if SEQAN_ENABLE_TESTING

// variadic macros are not supported by VS 2003 and before
#if !defined(_MSC_VER) || (_MSC_VER >= 1400)

#if SEQAN_ENABLE_DEBUG

// Force a test failure.
//
// Usage:  SEQAN_ASSERT_FAIL("Failed at position %d", pos);
#define SEQAN_ASSERT_FAIL(...)                                          \
    do {                                                                \
        ::seqan::ClassTest::forceFail(__FILE__, __LINE__,               \
                                      __VA_ARGS__);                     \
        ::seqan::ClassTest::fail();                                     \
    } while (false)


// Equality assertion without a comment.
//
// Usage:  SEQAN_ASSERT_EQ(4, 4);
#define SEQAN_ASSERT_EQ(_arg1, _arg2)                                   \
    do {                                                                \
        if (!::seqan::ClassTest::testEqual(__FILE__, __LINE__,          \
                                           (_arg1), #_arg1,             \
                                           (_arg2), #_arg2)) {          \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Equality assertion with a comment.
//
// Usage:  SEQAN_ASSERT_EQ(4, 4);
#define SEQAN_ASSERT_EQ_MSG(_arg1, _arg2, ...)                          \
    do {                                                                \
        if (!::seqan::ClassTest::testEqual(__FILE__, __LINE__,          \
                                           (_arg1), #_arg1,             \
                                           (_arg2), #_arg2,             \
                                           __VA_ARGS__)) {              \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// In-delta-environment assertion without a comment.
//
// Usage:  SEQAN_ASSERT_IN_DELTA(4.1, 4, 0.1);
#define SEQAN_ASSERT_IN_DELTA(_arg1, _arg2, _arg3)                      \
    do {                                                                \
        if (!::seqan::ClassTest::testInDelta(__FILE__, __LINE__,        \
                                             (_arg1), #_arg1,           \
                                             (_arg2), #_arg2,           \
                                             (_arg3), #_arg3)) {        \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// In-delta-environment assertion witha comment.
//
// Usage:  SEQAN_ASSERT_IN_DELTA_MSG(4.1, 4, 0.1, "3.9 <= 4.1 <= 4.1");
#define SEQAN_ASSERT_IN_DELTA_MSG(_arg1, _arg2, _arg3, ...)             \
    do {                                                                \
        if (!::seqan::ClassTest::testInDelta(__FILE__, __LINE__,        \
                                             (_arg1), #_arg1,           \
                                             (_arg2), #_arg2,           \
                                             (_arg3), #_arg3,           \
                                             __VA_ARGS__)) {            \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Inequality assertion without a comment.
//
// Usage:  SEQAN_ASSERT_NEQ(4, 5);
#define SEQAN_ASSERT_NEQ(_arg1, _arg2)                                  \
    do {                                                                \
        if (!::seqan::ClassTest::testNotEqual(__FILE__, __LINE__,       \
                                              (_arg1), #_arg1,          \
                                              (_arg2), #_arg2)) {       \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Inequality assertion with a comment.
//
// Usage:  SEQAN_ASSERT_NEQ(4, 5);
#define SEQAN_ASSERT_NEQ_MSG(_arg1, _arg2, ...)                         \
    do {                                                                \
        if (!::seqan::ClassTest::testNotEqual(__FILE__, __LINE__,       \
                                              (_arg1), #_arg1,          \
                                              (_arg2), #_arg2,          \
                                              __VA_ARGS__)) {           \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Less-than-or-equal assertion without a comment.
#define SEQAN_ASSERT_LEQ(_arg1, _arg2)                                  \
    do {                                                                \
        if (!::seqan::ClassTest::testLeq(__FILE__, __LINE__,            \
                                         (_arg1), #_arg1,               \
                                         (_arg2), #_arg2)) {            \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Less-than-or-equal assertion with a comment.
#define SEQAN_ASSERT_LEQ_MSG(_arg1, _arg2, ...)                         \
    do {                                                                \
        if (!::seqan::ClassTest::testLeq(__FILE__, __LINE__,            \
                                         (_arg1), #_arg1,               \
                                         (_arg2), #_arg2,               \
                                         __VA_ARGS__)) {                \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Less-than assertion without a comment.
#define SEQAN_ASSERT_LT(_arg1, _arg2)                                   \
    do {                                                                \
        if (!::seqan::ClassTest::testLt(__FILE__, __LINE__,             \
                                        (_arg1), #_arg1,                \
                                        (_arg2), #_arg2)) {             \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Less-than assertion with a comment.
#define SEQAN_ASSERT_LT_MSG(_arg1, _arg2, ...)                          \
    do {                                                                \
        if (!::seqan::ClassTest::testLt(__FILE__, __LINE__,             \
                                        (_arg1), #_arg1,                \
                                        (_arg2), #_arg2,                \
                                        __VA_ARGS__)) {                 \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Greater-than-or-equal assertion without a comment.
#define SEQAN_ASSERT_GEQ(_arg1, _arg2)                                  \
    do {                                                                \
        if (!::seqan::ClassTest::testGeq(__FILE__, __LINE__,            \
                                         (_arg1), #_arg1,               \
                                         (_arg2), #_arg2)) {            \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Greater-than-or-equal assertion with a comment.
#define SEQAN_ASSERT_GEQ_MSG(_arg1, _arg2, ...)                         \
    do {                                                                \
        if (!::seqan::ClassTest::testGeq(__FILE__, __LINE__,            \
                                         (_arg1), #_arg1,               \
                                         (_arg2), #_arg2,               \
                                         __VA_ARGS__)) {                \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Greater-than assertion without a comment.
#define SEQAN_ASSERT_GT(_arg1, _arg2)                                   \
    do {                                                                \
        if (!::seqan::ClassTest::testGt(__FILE__, __LINE__,             \
                                        (_arg1), #_arg1,                \
                                        (_arg2), #_arg2)) {             \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Greater-than assertion with a comment.
#define SEQAN_ASSERT_GT_MSG(_arg1, _arg2, ...)                          \
    do {                                                                \
        if (!::seqan::ClassTest::testGt(__FILE__, __LINE__,             \
                                        (_arg1), #_arg1,                \
                                        (_arg2), #_arg2,                \
                                        __VA_ARGS__)) {                 \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// TODO(holtgrew): Rename to SEQAN_ASSERT_TRUE once that name is free.;
// Trueness assertion with a comment.
//
// Usage:  SEQAN_ASSERT_TRUE(false);
#define SEQAN_ASSERT_TRUE(_arg1)                                        \
    do {                                                                \
        if (!::seqan::ClassTest::testTrue(__FILE__, __LINE__,           \
                                          (_arg1), #_arg1)) {           \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// TODO(holtgrew): Rename to SEQAN_ASSERT_TRUE once that name is free.;
// Trueness assertion with a comment.
#define SEQAN_ASSERT_TRUE_MSG(_arg1, ...)                               \
    do {                                                                \
        if (!::seqan::ClassTest::testTrue(__FILE__, __LINE__,           \
                                          (_arg1), #_arg1,              \
                                          __VA_ARGS__)) {             \
            ::seqan::ClassTest::fail();                                 \
        }                                                               \
    } while (false)


// Falseness assertion without a comment.
//
// Usage:  SEQAN_ASSERT_NOT(false);
#define SEQAN_ASSERT_NOT(_arg1)                                       \
    do {                                                              \
        if (!::seqan::ClassTest::testFalse(__FILE__, __LINE__,        \
                                           (_arg1), #_arg1)) {        \
            ::seqan::ClassTest::fail();                               \
        }                                                             \
    } while (false)


// Falseness assertion with a comment.
#define SEQAN_ASSERT_NOT_MSG(_arg1, ...)                              \
    do {                                                              \
        if (!::seqan::ClassTest::testFalse(__FILE__, __LINE__,        \
                                           (_arg1), #_arg1,           \
                                           __VA_ARGS__)) {          \
            ::seqan::ClassTest::fail();                               \
        }                                                             \
    } while (false)


#else  // #if SEQAN_ENABLE_DEBUG

#define SEQAN_ASSERT_EQ(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_EQ_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_NEQ(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_NEQ_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_LEQ(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_LEQ_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_LT(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_LT_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_GEQ(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_GEQ_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_GT(_arg1, _arg2) do {} while (false)
#define SEQAN_ASSERT_GT_MSG(_arg1, _arg2, ...) do {} while (false)
#define SEQAN_ASSERT_TRUE(_arg1) do {} while (false)
#define SEQAN_ASSERT_TRUE_MSG(_arg1, ...) do {} while (false)
#define SEQAN_ASSERT_NOT(_arg1) do {} while (false)
#define SEQAN_ASSERT_NOT_MSG(_arg1, ...) do {} while (false)
#define SEQAN_ASSERT_FAIL(...) do {} while (false)

#endif  // #if SEQAN_ENABLE_DEBUG

#else // no variadic macros

#if SEQAN_ENABLE_DEBUG
inline void SEQAN_ASSERT_FAIL(const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	::seqan::ClassTest::vforceFail("", 0, comment, args);
    ::seqan::ClassTest::fail();
	va_end(args);
}

template <typename T1, typename T2, typename T3>
void SEQAN_ASSERT_IN_DELTA(T1 const &_arg1, T2 const &_arg2, T3 const &_arg3)
{
	if (!::seqan::ClassTest::testInDelta("", 0, _arg1, "", _arg2, "", _arg3, ""))
		::seqan::ClassTest::fail();
}

template <typename T1, typename T2, typename T3>
void SEQAN_ASSERT_IN_DELTA_MSG(T1 const &_arg1, T2 const &_arg2, T3 const &_arg3, const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	if (!::seqan::ClassTest::vtestInDelta("", 0, _arg1, "", _arg2, "", _arg3, "", comment, args))
		::seqan::ClassTest::fail();
	va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_EQ(T1 const &_arg1, T2 const &_arg2)
{
	if (!::seqan::ClassTest::testEqual("", 0, _arg1, "", _arg2, ""))
		::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_EQ_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	if (!::seqan::ClassTest::vtestEqual("", 0, _arg1, "", _arg2, "", comment, args))
		::seqan::ClassTest::fail();
	va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_NEQ(T1 const &_arg1, T2 const &_arg2)
{
	if (!::seqan::ClassTest::testNotEqual("", _arg1, "", _arg2, ""))
		::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_NEQ_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	if (!::seqan::ClassTest::vtestNotEqual("", _arg1, "", _arg2, "", comment, args))
		::seqan::ClassTest::fail();
	va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_LEQ(T1 const &_arg1, T2 const &_arg2)
{
	if (!::seqan::ClassTest::testLeq("", 0, _arg1, "", _arg2, ""))
		::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_LEQ_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	if (!::seqan::ClassTest::vtestLeq("", 0, _arg1, "", _arg2, "", comment, args))
		::seqan::ClassTest::fail();
	va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_LT(T1 const &_arg1, T2 const &_arg2)
{
	if (!::seqan::ClassTest::testLt("", 0, _arg1, "", _arg2, ""))
		::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_LT_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	if (!::seqan::ClassTest::vtestLt("", 0, _arg1, "", _arg2, "", comment, args))
		::seqan::ClassTest::fail();
	va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_GEQ(T1 const &_arg1, T2 const &_arg2)
{
	if (!::seqan::ClassTest::testGeq("", 0, _arg1, "", _arg2, ""))
		::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_GEQ_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	if (!::seqan::ClassTest::vtestGeq("", 0, _arg1, "", _arg2, "", comment, args))
		::seqan::ClassTest::fail();
	va_end(args);
}

template <typename T1, typename T2>
void SEQAN_ASSERT_GT(T1 const &_arg1, T2 const &_arg2)
{
	if (!::seqan::ClassTest::testGt("", 0, _arg1, "", _arg2, ""))
		::seqan::ClassTest::fail();
}

template <typename T1, typename T2>
void SEQAN_ASSERT_GT_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	if (!::seqan::ClassTest::vtestGt("", 0, _arg1, "", _arg2, "", comment, args))
		::seqan::ClassTest::fail();
	va_end(args);
}

template <typename T1>
void SEQAN_ASSERT_TRUE(T1 const &_arg1)
{
	if (!::seqan::ClassTest::testTrue("", 0, _arg1, ""))
		::seqan::ClassTest::fail();
}

template <typename T1>
void SEQAN_ASSERT_TRUE_MSG(T1 const &_arg1, const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	if (!::seqan::ClassTest::vtestTrue("", 0, _arg1, "", comment, args))
		::seqan::ClassTest::fail();
	va_end(args);
}

template <typename T1>
void SEQAN_ASSERT_NOT(T1 const &_arg1)
{
	if (!::seqan::ClassTest::testFalse("", 0, _arg1, ""))
		::seqan::ClassTest::fail();
}

template <typename T1>
void SEQAN_ASSERT_NOT_MSG(T1 const &_arg1, const char *comment, ...)
{
	va_list args;
	va_start(args, comment);
	if (!::seqan::ClassTest::vtestFalse("", 0, _arg1, "", comment, args))
		::seqan::ClassTest::fail();
	va_end(args);
}

#else // #if SEQAN_ENABLE_DEBUG

inline void SEQAN_ASSERT_FAIL(const char *comment, ...) {}
template <typename T1, typename T2, typename T3> void SEQAN_ASSERT_IN_DELTA(T1 const &_arg1, T2 const &_arg2, T3 const &_arg3) {}
template <typename T1, typename T2, typename T3> void SEQAN_ASSERT_IN_DELTA_MSG(T1 const &_arg1, T2 const &_arg2, T3 const &_arg3, const char *comment, ...) {}
template <typename T1, typename T2> void SEQAN_ASSERT_EQ(T1 const &_arg1, T2 const &_arg2) {}
template <typename T1, typename T2> void SEQAN_ASSERT_EQ_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...) {}
template <typename T1, typename T2> void SEQAN_ASSERT_NEQ(T1 const &_arg1, T2 const &_arg2) {}
template <typename T1, typename T2> void SEQAN_ASSERT_NEQ_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...) {}
template <typename T1, typename T2> void SEQAN_ASSERT_LEQ(T1 const &_arg1, T2 const &_arg2) {}
template <typename T1, typename T2> void SEQAN_ASSERT_LEQ_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...) {}
template <typename T1, typename T2> void SEQAN_ASSERT_LT(T1 const &_arg1, T2 const &_arg2) {}
template <typename T1, typename T2> void SEQAN_ASSERT_LT_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...) {}
template <typename T1, typename T2> void SEQAN_ASSERT_GEQ(T1 const &_arg1, T2 const &_arg2) {}
template <typename T1, typename T2> void SEQAN_ASSERT_GEQ_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...) {}
template <typename T1, typename T2> void SEQAN_ASSERT_GT(T1 const &_arg1, T2 const &_arg2) {}
template <typename T1, typename T2> void SEQAN_ASSERT_GT_MSG(T1 const &_arg1, T2 const &_arg2, const char *comment, ...) {}
template <typename T1> void SEQAN_ASSERT_TRUE(T1 const &_arg1) {}
template <typename T1> void SEQAN_ASSERT_TRUE_MSG(T1 const &_arg1, const char *comment, ...) {}
template <typename T1> void SEQAN_ASSERT_NOT(T1 const &_arg1) {}
template <typename T1> void SEQAN_ASSERT_NOT_MSG(T1 const &_arg1, const char *comment, ...) {}

#endif // #if SEQAN_ENABLE_DEBUG

#endif // no variadic macros

// Returns a string (of type char*) with the path to the called binary.
//
// Use this to locate files relative to the test binary.
#define SEQAN_PROGRAM_PATH                      \
    ::seqan::ClassTest::StaticData::basePath()

// Returns a const char * string with the path to the projects directory.
#define SEQAN_PATH_TO_PROJECTS()                        \
    ::seqan::ClassTest::StaticData::pathToProjects()


// Returns the POSIX int file handle to an open file.
// TODO(holtgrewe): Uncomment if openTempFile has been implemented.
// #define SEQAN_OPEN_TEMP_FILE() (::seqan::ClassTest::openTempFile())


// Returns a temporary filename.
#define SEQAN_TEMP_FILENAME() (::seqan::ClassTest::tempFileName())


#if SEQAN_ENABLE_CHECKPOINTS

// Create a check point at the point where the macro is placed.
// TODO(holtgrew): Should be called SEQAN_CHECK_POINT to be consistent.
#define SEQAN_CHECKPOINT                                        \
    ::seqan::ClassTest::registerCheckPoint(__LINE__, __FILE__);

// Call the check point verification code for the given file.
#define SEQAN_VERIFY_CHECKPOINTS(filename)          \
    ::seqan::ClassTest::verifyCheckPoints(filename)

#else  // #if SEQAN_ENABLE_CHECKPOINTS

#define SEQAN_CHECKPOINT

// If checkpoints are to be verified if testing is disabled then print
// a warning.
#define SEQAN_VERIFY_CHECKPOINTS(filename)                              \
    do {                                                                \
        fprintf(stderr, ("WARNING: Check point verification is "        \
                         "disabled. Trying to verify %s from %s:%d.\n"), \
                filename, __FILE__, __LINE__);                          \
    } while(false)

#endif  // #if SEQAN_ENABLE_CHECKPOINTS

#if !SEQAN_ENABLE_TESTING

#define SEQAN_BEGIN_TESTSUITE(suite_name)                               \
    int main(int argc, char **argv) {                                   \
    (void) argc;                                                        \
    (void) argv;                                                        \
    fprintf(stderr, "Warning: SEQAN_ENABLE_TESTING is wrong and you used the macro SEQAN_BEGIN_TESTSUITE!\n");
#define SEQAN_END_TESTSUITE }
#define SEQAN_CALL_TEST(test_name) do { SEQAN_TEST_ ## test_name(); } while (false)
#define SEQAN_SKIP_TEST do {} while (false)

#endif  // #if !SEQAN_ENABLE_TESTING

}  // namespace seqan

#endif  // SEQAN_BASIC_BASIC_TESTING_H_
