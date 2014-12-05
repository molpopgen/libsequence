#Unit tests for libsequence

##Dependencies

1. Make sure that libsequence is compiled in the parent directory
2. The [boost](http://boost.org) unit testing library is used by these tests.  Currently, autoconf does _not_ check for this dependency.  Make sure that the library is installed

###A word of caution

I develop the library and the tests on an Ubuntu Linux machine.  The library is written in C++11 and tested primarily using GCC and secondarily using clang++ (the default compiler on current-era OS X).  On Ubuntu 14.04, I have observed that the unit tests fail to compile with clang++.  I have no tracked down if this due to the Ubuntu boost packages being compiled with GCC, without C++11 awareness, some issue with clang++ and boost's unit testing library, or some complex interaction amongst those possibilities.  However, I have confirmed that the unit testing compiles and works find on OS X Yosemite using clang++.

##Compiling the tests

```
make check
```

##Running the tests

```
sh runTests.sh
```

If you really want all the details, then execute this instead:

```
BOOST_TEST_LOG_LEVEL=all sh runTests.sh
```

The boost unit testing library will report any errors in any testing modules.

Note that some tests may intentionally cause errors.  When that it the case, a message stating that the error is intentional will appear on screen along with the error.

##Notes

* The tests are statically-linked against the version of libsequence compiled in the parent directory.  This is done so that there is no confusion that the tests are testing the current code, and not some other version of the library installed on your system.
* More tests are needed!