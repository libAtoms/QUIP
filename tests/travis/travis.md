Travis CI Build
============

QUIP makes use of [Travis CI](travisci.org) for running build testing and regression tests for pushed commits. All of the Travis specific files are in the `tests/travis` directory. Some of the files need to be edited to build all features successfully:

Setup
-------

The build should work as-is for basic builds. On push requests only `VANILLA` configurations will be built and only for the default compiler set, all other configurations will run but not do anything (so appearing as successful jobs). Daily cron jobs will build the `ALL` configurations (i.e. including all modules like GAP and ThirdParty) and refresh the documentation (setup for crons is on the Travis website). A weekly cron job builds everything for large range of gcc versions.

To set up the build to work for GAP and pushing the docs the [Travis Client](https://github.com/travis-ci/travis.rb) Is needed to encrypt sensitive data. Install this and work from inside a clone of the repository.

Some of the dependencies are hosted elsewhere, with limited access. It is assumed that these can be retrieved over SSH with an SSH key.

* Copy `get_deps.sh.src` to `get_deps.sh`.
* Add `src_url` so that it points to the where the dependencies are and put the contents of your private key in `src_key`, including the newlines: 
```bash
src_url="git@server.ac.uk"
src_key="-----BEGIN RSA PRIVATE KEY-----
Y253ZWlvbmlvZXduaW92cncgZXcgd2UgaDBmOHdlIG5pb2Z...
...
-----END RSA PRIVATE KEY-----"
```
* Encrypt the file: `travis encrypt-file get_deps.sh`
* Change the `openssl ...` line in `.travis.yml` if necessary.
* Add the `get_deps.sh.enc` changes to the repository; DO NOT ADD `get_deps.sh`!
* Commit and push to check if building `ALL` works

To build the documentation, Travis needs to be able to push to the gh-pages branch of the repository. Create an access token from https://github.com/settings/tokens to use (these can be easily revoked from the same page if something goes wrong). Replace username and token with someone who can push to the repository. These values will be masked in the build.

```console
$ travis encrypt "PAGES_URL=https://<username>:<token>@github.com/libAtoms/QUIP" 

```

Replace the `secure: "..."` line in the `.travis.yml` file with the new value. Commit and push to test if building the docs works. There must be a value of `env` with `DOCS=true` to build the docs.

Compilers
------------

The current default compiler on Travis is `gcc-4.6`, the weekly build  also tests `gcc-4.4` `gcc-4.9`, `gcc-5` and `gcc-6` packages from `ubuntu-toolchain-r-test`.  More compilers can be added by adding their `apt` packages and adding the appropriate `GCC_VERSION` value in `env`. Intel compiler compatibility is tested elsewhere.

Tests are also run against Python 2.7 only (2.6 support has been dropped Feb 2017).

Builds use the `sudo` enabled runners as installing each compiler is one of the slowest steps of the build. Using sudo allows only the required compilers to be installed for the specific build and to skip dependencies altogether for skipped builds.

Procedure
------------

Several steps are carried out in the build and testing:

* Check to see if the sub-job should run
* Install dependencies
* Set the compilers according to the value of `GCC_VERSION`
* Pull the dependencies, if required
* make
* make libquip
* install ase
* make quippy
* make test
* Build the docs (including atomeye), for one of the builds.

