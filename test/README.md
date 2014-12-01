#Unit tests for libsequence

##Dependencies

1. Make sure that libsequence is compiled in the parent directory
2. The [boost](http://boost.org) unit testing library is used by these tests.  Currently, autoconf does _not_ check for this dependency.  Make sure that the library is installed

##Compiling the tests

```
make check
```

##Running the tests

```
sh runTests.sh
```

The boost unit testing library will report any errors in any testing modules.

Note that some tests may intentionally cause errors.  When that it the case, a message stating that the error is intentional will appear on screen along with the error.

##Notes

* The tests are statically-linked against the version of libsequence compiled in the parent directory.  This is done so that there is no confusion that the tests are testing the current code, and not some other version of the library installed on your system.
* More tests are needed!