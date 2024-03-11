Releases
========

1.1.3
-----

Bug fix:
^^^^^^^^

- In 1d when kde could not be computing with recent version of scipy. baredSC was stoped raising a `linalg.LinAlgError`.

1.1.2
-----

Bug fix:
^^^^^^^^

- When ``--minNeff`` was set with ``--nsampMCMC``, the number of samples was not multiplied more than once by 10.

Improvements:
^^^^^^^^^^^^^

- The output ``means.txt.gz`` is now in documentation.
- Big linting effort.


1.1.1
-----

Improvements:
^^^^^^^^^^^^^

- The online documentation has been improved.

- More information obtained by ``--help``.


1.1.0
-----

Improvements:
^^^^^^^^^^^^^

- Support annData as input with ``--inputAnnData``


1.0.0
-----

First release
