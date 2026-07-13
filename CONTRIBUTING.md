# Contributing to VariantValidator

Thank you for your interest in contributing to VariantValidator.

VariantValidator is a community-driven open-source project. We welcome contributions that improve the software, documentation and user experience. Whether you are fixing a bug, implementing a new feature or improving the documentation, we encourage contributors to work with the community throughout the development process.

---

## Getting Started

Before contributing, install a local development version of VariantValidator by following the documentation.

See the **Installation Guide**:

- `docs/installation/index.md`

To ensure a consistent development environment across all contributors and continuous integration (CI) systems, **all development must be carried out using the supplied Conda environment**. Contributions developed outside the supported environment may not be accepted if they introduce compatibility or reproducibility issues.

---

## Development Workflow

All development follows the Git branching workflow.

1. Fork the `openvar/variantValidator` repository.
2. Create a feature branch from the **develop** branch.
3. Make your changes.
4. Add or update tests.
5. Ensure the complete test suite passes.
6. Push your branch to your fork.
7. Open a Pull Request from your fork into the **develop** branch of the `openvar/variantValidator` repository.

> **Important**
>
> - Always create feature branches from **develop**.
> - Always submit Pull Requests into **develop**.
> - Pull Requests targeting **master** will not be accepted.

---

## Planning New Features

VariantValidator is developed in collaboration with its user community.

Before implementing a new feature, please open a GitHub Issue describing the proposed functionality. This allows discussion of the proposed design before development begins and helps ensure that new functionality meets community requirements.

Bug fixes may be submitted directly if they are straightforward, although opening an Issue beforehand is still encouraged.

---

## Testing Requirements

All contributions must include appropriate automated tests.

### New Functionality

New functionality should be accompanied by automated tests providing **at least 90% code coverage** of the new or modified code.

### Bug Fixes

Every bug fix must include a **regression test** demonstrating the original problem and verifying that it remains fixed.

Where practical, VariantValidator strongly prefers functional and integration tests over mocked tests, as these more accurately reflect real-world behaviour.

---

## Running the Test Suite

Before submitting a Pull Request, ensure the complete test suite passes.

If using the Docker database installation:

```bash
pytest -n 4
```

If using a fully native database installation:

```bash
pytest
```

All tests must pass before a Pull Request is submitted.

---

## Coding Standards

Please aim to:

- Follow the existing coding style.
- Write clear, readable and well-documented code.
- Prefer small, self-contained functions with clear docstrings.
- Keep Pull Requests focused on a single logical change.
- Add or update documentation where appropriate.
- Avoid introducing unnecessary dependencies.

---

## Formatting

Please ensure all contributed code follows the existing project style.

In particular:

- Wrap source code at approximately **79–88 characters** where practical.
- Wrap docstrings and comments at approximately **72 characters** where practical.
- Use four spaces for indentation (never tabs).
- Use UTF-8 encoded source files.
- Follow existing naming conventions and project structure.
- Avoid unnecessary whitespace-only changes.
- Run formatting tools where appropriate before committing.

Where an existing file follows a particular style, please maintain that style for consistency rather than introducing a different formatting convention.

---

## Commit Messages

Please write clear, descriptive commit messages.

The first line (the commit summary) should:

- Be written in the imperative mood (for example, *Add*, *Fix*, *Update* or *Refactor*).
- Be no longer than **50 characters** where practical.
- Summarise the change without ending in a full stop.

If additional explanation is required, leave a blank line after the summary and provide a more detailed description wrapped at approximately **72 characters**.

For example:

```text
Add hgvs2reference Python API

Document the hgvs2reference Python API, including
supported inputs, return values and examples.
```

Small, logically related commits are preferred over large commits containing unrelated changes, as they simplify code review and future maintenance.

---

## Submitting a Pull Request

When your work is complete:

1. Push your feature branch to your fork.
2. Open a Pull Request into the **develop** branch of the `openvar/variantValidator` repository.
3. Reference any related GitHub Issues.
4. Provide a clear description of the changes.
5. Ensure all Continuous Integration (CI) checks complete successfully.

Where practical, each Pull Request should address a single feature, enhancement or bug fix.

All Pull Requests are reviewed by the VariantValidator maintainers before being merged into **develop**. Changes are promoted to **master** as part of the normal release process.

---

## Reporting Bugs

Bug reports are always welcome.

Please include:

- VariantValidator version
- Operating system
- Installation method
- Example input
- Expected behaviour
- Observed behaviour
- Any relevant error messages or traceback

Providing a minimal reproducible example greatly helps us investigate issues.

---

## Documentation

Documentation is considered part of the software.

New features and behavioural changes should include corresponding documentation updates where appropriate.

Documentation improvements are always appreciated. If you notice unclear explanations, missing examples or broken links, please feel free to submit a Pull Request.

---

Thank you for helping improve VariantValidator.
