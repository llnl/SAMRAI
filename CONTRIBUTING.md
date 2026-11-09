# Contributing to SAMRAI

This document is intended for developers who want to add new features or
bug fixes to SAMRAI. It assumes you have some familiarity with Git and GitHub. It
explains what a good pull request (PR) looks like, and the tests your
PR must pass before it can be merged into SAMRAI.

## Forking SAMRAI

If you aren't a SAMRAI developer at LLNL, you won't have permission to push
new branches to the repository. First, you should create a
[fork](https://github.com/LLNL/SAMRAI#fork-destination-box). This will create a
copy of the SAMRAI repository that you own, and will ensure you can push your
changes to GitHub and create pull requests.

## Developing a New Feature

New features should be based on the `master` branch. When you want to create a
new feature, first ensure you have an up-to-date copy of the `master` branch:

    $ git checkout master
    $ git pull origin master

You can now create a new branch to develop your feature on:

    $ git checkout -b feature/<name-of-feature>

Proceed to develop your feature on this branch, and add tests that will exercise
your new code. If you are creating new methods or classes, please add Doxygen
documentation.

Once your feature is complete and your tests are passing, you can push your
branch to GitHub and create a PR.

## Developing a Bug Fix

First, check if the change you want to make has already been fixed in `master`. If so,
we suggest you either start using the `master` branch or temporarily apply the
fix to whichever version of SAMRAI you are using.

If the bug remains unsolved, make sure you have an up-to-date copy
of the `develop` branch:

    $ git checkout develop
    $ git pull origin develop

Then create a new branch for your bug fix:

    $ git checkout -b bugfix/<name-of-bug>

Start by adding a test that reproduces the bug you have found. Then, develop your
bug fix as normal, and ensure to run `make test` to check that your changes actually fix
the bug.

Once you're finished, push your branch to GitHub and create a PR.

## Creating a Pull Request

You can create a new PR [here](https://github.com/LLNL/SAMRAI/compare). GitHub
has a good [guide](https://help.github.com/articles/about-pull-requests/) to PR
basics if you need more information. Ensure that your PR base is the
`develop` branch of SAMRAI.

Add a descriptive title explaining the bug you fixed or the feature you have
added, and provide a detailed description of the changes you have made in the comment
box.

Once your PR is created, it will be run through our automated tests and
also reviewed by SAMRAI team members. If the branch passes both the
tests and the review, it will be merged into SAMRAI.

## Tests

SAMRAI uses Travis CI for continuous integration tests. Our tests are
automatically run against every new pull request, and passing all tests is a
requirement for merging your PR. If you are developing a bug fix or a new
feature, please add a test to verify the correctness of your new code. SAMRAI
is used on a wide variety of systems with different configurations, and adding
new tests helps ensure that all features work as expected across these
environments.

SAMRAI's tests are all in the `/source/test` directory and are split up by module to match different parts of the codebase. Please review existing tests to determine where your test should go.
