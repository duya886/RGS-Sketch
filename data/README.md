## **Notification**

**The datasets in our experiments are [CAIDA-2016](http://www.caida.org/data/passive/passive_2016_dataset.xml), [CAIDA-2019](http://www.caida.org/data/passive/passive_2019_dataset.xml), [eCommerce behavior data from a multi category store](https://www.kaggle.com/datasets/mkechinov/ecommerce-behavior-data-from-multi-category-store), and [eCommerce events history in a cosmetics shop](https://www.kaggle.com/datasets/mkechinov/ecommerce-events-history-in-cosmetics-shop).** We do not provide CAIDA datasets due to the license limit. We provide the binary data files of the two eCommerce datasets. Besides, we provide a Python code file to show how we encode the <flowid, element> records into our input binary data files.

## **Data Format**

The data files we use are binary. For a data item where the flowId and element are 32-bit integers, the data item is organized into a 64-bit integer: (element<<32)+flowId. Then, we write the data item into the binary file in little-endian.

For example, if the flowid is 1, and the element id is 3. The data item will be organized into 0x0000000300000001, and the value written into the binary data file is b'\x01\x00\x00\x00\x03\x00\x00\x00'. Then, when reading 8 bytes once in a C++ program (the byte ordering is little-endian), we will get 0x0000000300000001, where the higher part is the element, and the lower part is the flowId.

