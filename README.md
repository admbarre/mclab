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
- Create general interface to website
  - Make ordering/downloading into selectable options from this interface
- Batch updates to google sheets to avoid rate limiting issues
- Create GUI interface
  - Allow for specifying different users
  - Allow for clearer ordering
- Revamp download code to bring in everything not just latest and organize by date (this is important for reruns)
