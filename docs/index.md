# VariantValidator Documentation
[![VariantValidator CI](https://github.com/openvar/variantValidator/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/openvar/variantValidator/actions/workflows/ci.yml) [![codecov](https://codecov.io/github/openvar/variantValidator/graph/badge.svg?token=QWTxw5kiY4)](https://codecov.io/github/openvar/variantValidator)

Welcome to the VariantValidator documentation.

VariantValidator is an open-source software suite for validating, mapping, normalising and formatting genetic sequence variant descriptions according to the recommendations of the Human Genome Variation Society (HGVS).

The software is designed for researchers, clinical scientists, diagnostic laboratories, bioinformaticians and software developers working with genomic variation.

In addition to validating standard HGVS sequence variant descriptions, VariantValidator recognises many commonly encountered non-HGVS formats, legacy representations and common formatting mistakes. Where possible, these are automatically converted into valid HGVS sequence variant descriptions with informative warnings describing any corrections made.

The VariantValidator software suite also provides tools for variant formatting, transcript discovery and selection, and reference sequence retrieval.

---

## Documentation

The documentation is organised into two main sections.

### Installation

Install, configure and verify VariantValidator.

Topics include:

- Linux and macOS installation
- Windows installation using WSL2
- Docker Quick Start
- Full Docker deployment
- Configuration
- Installation verification
- Troubleshooting

→ [Installation Guide](installation/index.md)

---

### User Manual

Learn how to use VariantValidator and the accompanying tools through the command-line interface and Python API.

Topics include:

- VariantValidator
- VariantFormatter
- gene2transcripts
- hgvs2reference
- Supported input formats
- Output formats
- Transcript selection
- Python API
- Command-line tools

→ [User Manual](user-manual/index.md)

---

## Community

VariantValidator is a community-driven open-source project.

Bug reports, feature requests and suggestions for improvements are welcomed. Many new features originate directly from feedback from researchers, clinicians and diagnostic laboratories using VariantValidator in research and clinical practice.

If you encounter a problem or have a suggestion:

- Search the existing GitHub Issues.
- Open a new GitHub Issue to report a bug or request a feature.
- Contact the development team via the VariantValidator support page.

---

---

## Terms and conditions of use

### By continuing to use VariantValidator, you accept the following terms:
All contents of VariantValidator are protected by local and international copyright laws. User inputs are submitted for the purpose of interpreting genetic and clinical information (including clinical genomic and genetic data). Outputs returned may or may not have a causal association with disease phenotypes, irrespective of stated classifications or other information presented by VariantValidator. All information returned by VariantValidator, including variant classifications, is subject to change and there is no warranty, express or implied, as to its accuracy, completeness, or fitness for a particular purpose. Use of VariantValidator and information is subject to User responsibility and discretion. Clinical decisions regarding individual patient care should be carried out in conjunction with a healthcare professional with expertise in the relevant genes and diseases. We do not accept any liability for any injury, loss or damage incurred by use of, or reliance on, the information provided by VariantValidator.

## Web-services additional Terms and Conditions of use

---

### By continuing to use VariantValidator web-services, you accept the following terms:
1. Secure/Remote Access. All access and use of VariantValidator must be made via a secure network.

2. Variations in Content. The University of Leicester and the University of Manchester (The Universities) reserve the right, in their reasonable and good faith discretion, to remove or modify materials accessed by VariantValidator or outputs provided by VariantValidator because such materials contain errors or could be subject to an infringement or other adverse claim by a third party.

3. Remedial Action. Without limiting the above, The Universities may suspend delivery of the Service if it can be reasonably shown that a User’s failure to comply with this Agreement may cause irreparable harm to them.

4. Service Level. The Universities will use reasonable efforts to provide access to the Service on a continuous 24/7 basis (except for regularly scheduled maintenance when Service may be suspended) and free from viruses or other harmful software. The Universities shall not be liable for any failure or delay or interruption in the Service or failure of any equipment or telecommunications resulting from any cause beyond the Universities reasonable control. User is responsible for providing all required information.

4. No Warranty. The Universities make no warranty that Service is error free or that the use thereof will be uninterrupted and User acknowledges and agrees that the existence of such errors shall not constitute a breach of this Agreement. The Universities disclaim all other warranties with respect to Service, either express or implied, including but not limited to any implied warranties relating to quality, fitness for any particular purpose or ability to achieve a particular result. Limitation of Liability. Universities total liability for any claims, losses, damages or expenses whatsoever and howsoever caused (even if caused by the Universities negligence and/or breach of contract) shall be limited to the price paid to the Universities for the products or services that are the subject of the User’s claim.

### Changes to our terms of service

If we find it necessary to change these terms and conditions, including our Privacy Policy and our Terms and Conditions of use, for example to add provisions for a new service, or to curtail problematic use cases, we will provide you with at least either 1 calendar month, or 4 weeks, of notice, before these changes become effective. This notice will be provided via your registered email address. You are responsible for keeping our records of your email address up-to-date, and checking your emails from us, in order to receive such notice in a timely manner. If you continue to use any of the services covered by our terms of service after the notice period, this will be considered to be consent to the new terms of service. This notice period only applies to the terms of service, and does not constrain our right to change the functions of the services we provide, particularly in the case of upgrades, corrections of problematic user behaviour, and other technical issues.

## Privacy Policy

---

### Why we ask you to create an account
VariantValidator is an Open Source project exclusively funded by grant income and generous contributions from the University of Manchester (UK) and the University of Leicester (UK). In order to apply for external funding and to retain funding from our host institutions we must demonstrate the wider impact of our software throughout the global genomics community. By providing us with some basic information, you are helping us to demonstrate the need for the software and this will help us to build a convincing case which will enable the continued maintenance and future development of this community driven resource.

### Your data
We will not share any personal data with external organisations. We may use data that you have submitted to build a profile of our users. However, any such profiles will not include any identifiable data _e.g._ user names or email addresses.

### Why do I need to provide an email address
We require an email address for several reasons. The primary reason is to allow us to validate your account, _i.e._ to ensure that the requested account has originated from a genuine email address. Our batch services require an email address so that results data can be returned to you. Therefore, by creating and verifying your account, we can ensure that your data will only ever be returned to you.

We will only contact you directly for one of five reasons:

1. If we receive a direct query from you
2. Via automated responses generated by our batch tools
3. To inform you about changes to our Terms and Conditions
4. To inform you of major changes to our services that may impact your existing workflows
5. If our error logs indicate that you might need help in submitting validation requests that will succeed

For any further information, or to ask for help or guidance, please contact us at admin@variantvalidator.org or via our [web form](https://variantvalidator.org/help/contact/).

---

## Cite VariantValidator

If VariantValidator contributes to published work, please cite the following publications.

Freeman PJ, Hart RK, Gretton LJ, Brookes AJ, Dalgleish R. **VariantValidator: Accurate validation, mapping and formatting of sequence variation descriptions.** *Human Mutation.* 2018;39:61–68. https://doi.org/10.1002/humu.23348

Freeman PJ, Wagstaff JF, Fokkema IFAC, *et al.* **Standardizing variant naming in literature with VariantValidator to increase diagnostic rates.** *Nature Genetics.* 2024;56:2284–2286. https://doi.org/10.1038/s41588-024-01938-w

---

## License

VariantValidator is released under the terms of the GNU Affero General Public License version 3 (AGPL-3.0).

See [LICENSE.txt](https://github.com/openvar/variantValidator/blob/master/LICENSE.txt) for the full license text.

> Copyright (C) 2016–2026 VariantValidator Contributors

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="https://github.com/i3hsInnovation/resources/blob/master/images/UoM_logo.jpg?raw=true" width="40%" align="left"/>
<img src="https://github.com/i3hsInnovation/resources/blob/master/images/UoL-Logo-Full-Colour.png?raw=true" width="40%" align="right" />
<br clear="both"/>