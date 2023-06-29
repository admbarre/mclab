# mclab
Tools for interfacing with McLab's Sequencing service

## Overview
### McLab
- Module for interfacing with McLab's excel order form

### McOrder
- Module for placing order using google sheets api

### McDown
- Module for donwloading data from McLab

### process_results
- Module for scanning sequence reads for barcodes

### barcode
- Module for scanning for scanning for barcodes

## Future
- Clean up barcode code
    - ISSUE: need to address key errors that occur if a barcode isn't in
      dictionary
- Create general interface to website
    - FEATURE: Make ordering/downloading into selectable options from this interface
    - FEATURE: Check online data files against local data files and download
    accordingly
    - FEATURE (sheets): create a sheet for sequencing results because we might
      have re-runs of sequencing
- Batch updates to google sheets to avoid rate limiting issues
- Create GUI interface
  - Allow for specifying different users
  - Allow for clearer ordering
- Allow for full ordering sequence (entry into McLab) from program
  - be really careful with this, actual money LOL
