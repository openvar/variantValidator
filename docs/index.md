<img src="./static/img/logos/VV_logo.png" width="20%" />

# VariantValidator Documentation

[![VariantValidator CI](https://github.com/openvar/variantValidator/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/openvar/variantValidator/actions/workflows/ci.yml)
[![codecov](https://codecov.io/github/openvar/variantValidator/graph/badge.svg?token=QWTxw5kiY4)](https://codecov.io/github/openvar/variantValidator)

Welcome to the VariantValidator documentation.

VariantValidator is an open-source software suite for validating, mapping, normalising and formatting genetic sequence variant descriptions according to the recommendations of the Human Genome Variation Society (HGVS).

The software is designed for researchers, clinical scientists, diagnostic laboratories, bioinformaticians and software developers working with genomic variation.

In addition to validating standard HGVS sequence variant descriptions, VariantValidator recognises many commonly encountered non-HGVS formats, legacy representations and common formatting mistakes. Where possible, these are automatically converted into valid HGVS sequence variant descriptions with informative warnings describing any corrections made.

The VariantValidator software suite also provides tools for variant formatting, transcript discovery and selection, and reference sequence retrieval.

---

## Getting started

The VariantValidator documentation is organised according to how you intend to use VariantValidator. Choose the guide that best matches your workflow.

### I'm using the VariantValidator website

If you use VariantValidator through your web browser, the **Web Interface Guide** is the best place to start. It is intended for clinicians, researchers, laboratory scientists and other users who do not need to install VariantValidator or write Python code.

**Most users should begin with this guide.**

Topics include:

- Variant Validator
- Batch Validator
- Gene2Transcripts
- HGVS2Reference
- LOVD HGVS Syntax Checker
- Creating an account
- Downloading and interpreting results

→ [Web Interface Guide](vvweb/index.md)

---

### I'm installing VariantValidator

If you plan to install VariantValidator locally or deploy your own instance, see the **Installation Guide**.

Topics include:

- Linux and macOS installation
- Windows installation using WSL2
- Docker Quick Start
- Full Docker deployment
- Configuration
- Installation verification
- Installing the LOVD HGVS Syntax Checker
- Troubleshooting

→ [Installation Guide](installation/index.md)

---

### I'm developing with VariantValidator

If you are using VariantValidator in scripts, command-line workflows, pipelines or software development, see the **User Manual**.

Topics include:

- VariantValidator CLI and Python API
- VariantFormatter CLI and Python API
- gene2transcripts CLI and Python API
- hgvs2reference CLI and Python API
- LOVD HGVS Syntax Checker CLI and Python API
- Supported input formats
- Output formats
- Transcript selection
- Errors and error codes

→ [User Manual](user-manual/index.md)

---

### I'm integrating VariantValidator into my own software

If you are developing applications, web services or automated workflows that communicate with VariantValidator programmatically, see the **REST API Guide**.

Topics include:

- REST VariantValidator
- REST VariantValidator2
- LOVD HGVS Syntax Checker REST API
- Authentication
- Available endpoints
- API documentation
- Integration guidance

→ [REST API Guide](rest-vv/index.md)

---

## Community Driven

<img src="./static/img/Our_Community.png" align="left"/>

*Image by [Rosaria](https://www.instagram.com/2drosaria/?hl=en-gb).*

VariantValidator is a community-driven, open-source project developed in close collaboration with researchers, clinicians, diagnostic laboratories, educators and software developers across the international genomics community.

We actively collaborate with complementary community resources to improve the interpretation, standardisation and exchange of genomic variant data. In particular, we work closely with the **LOVD (Leiden Open Variation Database)** development team, jointly developing and integrating software components that improve HGVS nomenclature support and the interoperability of genomic variant interpretation tools. This collaborative approach helps reduce duplicated effort, encourages shared standards, and benefits the wider genomics community.

We follow an agile development approach, releasing new features regularly and continually improving the platform based on community feedback. Many of the features available in VariantValidator today originated from suggestions made by users, while bug reports, feature requests and collaborations with partner projects help us prioritise future development.

If you discover a bug, have an idea for a new feature, or simply think something could work better, we'd love to hear from you. You can:

- Search existing issues on GitHub.
- Open a new GitHub issue to report a bug or request a feature.
- Contact the development team through our support page.

Community engagement and collaboration with partner projects play a vital role in the continued development of VariantValidator, helping us ensure the software keeps pace with the rapidly evolving field of genomic medicine.

---

## Terms and Conditions of Use

### By continuing to use VariantValidator, you accept the following terms:

All contents of VariantValidator are protected by local and international copyright laws. User inputs are submitted for the purpose of interpreting genetic and clinical information (including clinical genomic and genetic data). Outputs returned may or may not have a causal association with disease phenotypes, irrespective of stated classifications or other information presented by VariantValidator. All information returned by VariantValidator, including variant classifications, is subject to change and there is no warranty, express or implied, as to its accuracy, completeness, or fitness for a particular purpose. Use of VariantValidator and information is subject to user responsibility and discretion. Clinical decisions regarding individual patient care should be carried out in conjunction with a healthcare professional with expertise in the relevant genes and diseases. We do not accept any liability for any injury, loss or damage incurred by use of, or reliance on, the information provided by VariantValidator.

## Web-services Additional Terms and Conditions of Use

---

### By continuing to use VariantValidator web-services, you accept the following terms:

1. Secure/Remote Access. All access and use of VariantValidator must be made via a secure network.

2. Variations in Content. The University of Leicester and the University of Manchester (the Universities) reserve the right, in their reasonable and good faith discretion, to remove or modify materials accessed by VariantValidator or outputs provided by VariantValidator because such materials contain errors or could be subject to an infringement or other adverse claim by a third party.

3. Remedial Action. Without limiting the above, the Universities may suspend delivery of the Service if it can be reasonably shown that a User’s failure to comply with this Agreement may cause irreparable harm to them.

4. Service Level. The Universities will use reasonable efforts to provide access to the Service on a continuous 24/7 basis (except for regularly scheduled maintenance when the Service may be suspended) and free from viruses or other harmful software. The Universities shall not be liable for any failure, delay or interruption in the Service or failure of any equipment or telecommunications resulting from any cause beyond the Universities' reasonable control. The User is responsible for providing all required information.

5. No Warranty. The Universities make no warranty that the Service is error free or that its use will be uninterrupted. The User acknowledges and agrees that the existence of such errors shall not constitute a breach of this Agreement. The Universities disclaim all other warranties with respect to the Service, either express or implied, including but not limited to any implied warranties relating to quality, fitness for a particular purpose or ability to achieve a particular result.

6. Limitation of Liability. The Universities' total liability for any claims, losses, damages or expenses whatsoever and howsoever caused (including negligence and/or breach of contract) shall be limited to the price paid to the Universities for the products or services that are the subject of the User's claim.

### Changes to our Terms of Service

If we find it necessary to change these terms and conditions, including our Privacy Policy and our Terms and Conditions of Use, for example to add provisions for a new service or to curtail problematic use cases, we will provide you with at least one calendar month (or four weeks) notice before these changes become effective. This notice will be provided via your registered email address.

If you continue to use any of the services covered by our Terms of Service after the notice period, this will be considered consent to the updated terms.

---

## Privacy Policy

### Why we ask you to create an account

VariantValidator is an open-source project funded through grant income and contributions from the University of Manchester and the University of Leicester. Information provided during registration helps us demonstrate the impact of the software and supports future funding applications.

### Your data

We will not share personal data with external organisations. Aggregated, non-identifiable usage information may be used to understand how the software is being used.

### Why do I need to provide an email address?

We require an email address to:

1. Verify your account.
2. Return results generated by batch services.
3. Notify you of changes to our Terms and Conditions.
4. Notify you of major service changes that may affect existing workflows.
5. Contact you if our error logs indicate that we may be able to help resolve problems with submitted validation requests.

For further information or support, please contact us at **admin@variantvalidator.org** or via our [contact form](https://variantvalidator.org/help/contact/).

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

<img src="./static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="./static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>
