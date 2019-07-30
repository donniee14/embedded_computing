/*
    Donald Elmore and Kyle McMindes

    Huffman Codec
    
    References: GeeksForGeeks Huffman coding page was used for Huffman helper
        functions as they are widely known and used so no need to reinvent the wheel
*/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include "huffman_codec.h"

#define MAX_TREE_HT 256

unsigned int huffCodes[256];
int huffLengths[256] = {0};

//huffman tree node 
struct node {   
    unsigned char data; 
    unsigned freq; 
    struct node *left, *right; 
}; 
  
//collection of tree nodes
struct tree { 
    unsigned size; 
    unsigned capacity; 
    struct node** array; 
}; 
  
//create a new node
struct node* newNode(unsigned char data, unsigned freq) { 
    struct node* temp = (struct node*)malloc(sizeof(struct node)); 
    temp->left = temp->right = NULL; 
    temp->data = data; 
    temp->freq = freq; 
    return temp; 
}
  
//create a tree
struct tree* createTree(unsigned capacity) { 
    struct tree* tree = (struct tree*)malloc(sizeof(struct tree)); 
    tree->size = 0; 
    tree->capacity = capacity;  
    tree->array = (struct node**)malloc(tree->capacity * sizeof(struct node*)); 
    return tree; 
} 

//swap two nodes
void swapNodes(struct node** node1, struct node** node2) { 
    struct node* t = *node1; 
    *node1 = *node2; 
    *node2 = t; 
} 
  
//standard minheapify fn
void minHeapify(struct tree* tree, int idx) { 
    int smallest = idx; 
    int left = 2 * idx + 1; 
    int right = 2 * idx + 2; 
  
    if (left < tree->size && tree->array[left]->freq < tree->array[smallest]->freq) 
        smallest = left; 
  
    if (right < tree->size && tree->array[right]->freq < tree->array[smallest]->freq) 
        smallest = right; 
  
    if (smallest != idx) { 
        swapNodes(&tree->array[smallest], &tree->array[idx]); 
        minHeapify(tree, smallest); 
    } 
} 

//check is size one... duh
int isSizeOne(struct tree* T) { 
    return (T->size == 1); 
}
  
//find min from tree
struct node* findMin(struct tree* T) {
    struct node* temp = T->array[0]; 
    T->array[0] = T->array[T->size - 1]; 
    --T->size; 
    minHeapify(T, 0);
    return temp; 
} 
  
//insert a new node to tree
void insertTree(struct tree* tree, struct node* node) { 
    ++tree->size; 
    int i = tree->size - 1; 
    while (i && node->freq < tree->array[(i - 1) / 2]->freq) { 
        tree->array[i] = tree->array[(i - 1) / 2]; 
        i = (i - 1) / 2; 
    }
    tree->array[i] = node; 
} 
  
//build min heap
void buildTree(struct tree* T) { 
    int n = T->size - 1; 
    int i;  
    for (i = (n - 1) / 2; i >= 0; --i) 
        minHeapify(T, i); 
} 
  
//is the given node a leaf
int isLeaf(struct node* node) {  
    return !(node->left) && !(node->right);
} 

struct tree* makeTree(unsigned char data[], int freq[], int size) { 
    struct tree* tree = createTree(size); 

    for (int i = 0; i < size; ++i) 
        tree->array[i] = newNode(data[i], freq[i]); 
  
    tree->size = size; 
    buildTree(tree); 
    return tree; 
} 
  
//main function that builds Huffman tree
struct node* buildHuffmanTree(unsigned char data[], int freq[], int size) { 
    struct node *left, *right, *top; 
    struct tree* tree = makeTree(data, freq, size); 
  
    while (!isSizeOne(tree)) { 
        left = findMin(tree); 
        right = findMin(tree); 
        top = newNode('$', left->freq + right->freq); 
        top->left = left; 
        top->right = right; 
        insertTree(tree, top); 
    }
  
    return findMin(tree); 
} 
  
//print array of size n and save huffman codes and lengths
void printArr(unsigned char arr[], int n, struct node* root) { 
    int i;
    unsigned int byte = 0;
    for (i = 0; i < n; ++i) {
        //printf("%u", arr[i]);
        byte = byte | (arr[i] << (n-i-1));
    }
    huffCodes[(int)root->data] = byte;  //save code to codes array
    huffLengths[(int)root->data] = n;   //save code length to lengths array
    //printf("\n");
}

//print huffman codes
void outputCodes(struct node* root, unsigned char arr[], int top, int freq[]) { 
    if (root->left) { 
        arr[top] = 0; 
        outputCodes(root->left, arr, top + 1, freq); 
    }
    if (root->right) { 
        arr[top] = 1; 
        outputCodes(root->right, arr, top + 1, freq); 
    } 
    if (isLeaf(root) && (freq[(int)root->data] != 0)) {
        //printf("%u\t\t", root->data);
        printArr(arr, top, root);
    }
} 
  
//builds tree and traverses to print codes
void huff(unsigned char data[], int freq[], int size) { 
    struct node* root = buildHuffmanTree(data, freq, size); 
    unsigned char arr[MAX_TREE_HT] = {0}, top = 0;
    //printf("SYMBOL\t\tCODE\n");
    outputCodes(root, arr, top, freq); 
}

//compresses file
void compress(char *cmdArg) {
    printf("Compressing ...\n");
    FILE *fptin, *fptout;
    int *freqptr, i;

    //reading size of file
    fptin = fopen(cmdArg, "rb+");
    if(fptin == NULL)
        printf("Error opening file!\n");

    fseek(fptin, 0, SEEK_END);
    long int fileSize = ftell(fptin);
    fclose(fptin);

    //reading data to array of unsigned chars
    fptin = fopen(cmdArg, "rb+");
    if(fptin == NULL)
        printf("Error opening file!\n");

    unsigned char *in = (unsigned char *)malloc(fileSize);
    fread(in, sizeof(unsigned char), fileSize, fptin);
    fclose(fptin);

    fptin = fopen(cmdArg,"rb");

    if(fptin == NULL)
        printf("Error opening files!\n");
    
    //get the frequency of each character of the given file
    freqptr = getfreq(fptin);

    //initialize symbols array with byte values 0->255
    unsigned char byte, symbols[255];
    for (i = 0; i < 256; i++) {
        byte = i;
        symbols[i] = byte;
    }
    fclose(fptin);

    int size = sizeof(symbols) / sizeof(symbols[0]);
    huff(symbols, freqptr, size);

    /*
    printf("\nFINAL TABLE\n");
    printf("Symb\tFreq\tLen\tCode\n");
    for (i = 0; i < 256; i++) {
        printf("%d\t%d\t%d\t%u\n", symbols[i], freqptr[i], huffLengths[i], huffCodes[i]);
    } */

    //print dictionary
    //freq SPACE length SPACE code NEWLINE
    unsigned char space = ' ';
    unsigned char NL = '\n';
    fptout = fopen("dictionary", "wb");
    for (i = 0; i < 256; i++) {
        fwrite(&freqptr[i], sizeof(int), 1, fptout);
        fwrite(&space, 1, 1, fptout);
        fwrite(&huffLengths[i], sizeof(int), 1, fptout);
        fwrite(&space, 1, 1, fptout);
        fwrite(&huffCodes[i], sizeof(unsigned int), 1, fptout);
        fwrite(&NL, 1, 1, fptout);
    }
    fclose(fptout);

    //output compressed file
    unsigned char curr;
    int bufLength = 0;  //number of bits shifted into buf so far
    unsigned int buf1 = 0;
    fptout = fopen("compressed.huf", "wb");
    fptin = fopen(cmdArg, "rb");
    while (fread(&curr, 1, 1, fptin) == 1) {
        for (i = 0; i < 256; i++) {
            if (curr == symbols[i]) {
                //if space left in buf1, okay to write to buf1
                if ((32 - bufLength) >= huffLengths[i]) {    
                    buf1 = buf1 | (huffCodes[i] << (32 - bufLength - huffLengths[i]));
                    bufLength += huffLengths[i];
                    //if buffer fills perfectly
                    if (bufLength == 32) {
                        fwrite(&buf1, sizeof(unsigned int), 1, fptout);
                        bufLength = 0;
                        buf1 = 0;
                    }
                } else {
                    //not enough space in buf1
                    buf1 = buf1 | (huffCodes[i] >> (huffLengths[i] - (32 - bufLength)));
                    fwrite(&buf1, sizeof(unsigned int), 1, fptout);
                    buf1 = 0;
                    buf1 = buf1 | huffCodes[i] << (32 - (huffLengths[i] - (32 - bufLength)));
                    bufLength = (huffLengths[i] - (32 - bufLength));
                }
            }
        }
    }
    if (bufLength != 32) {
        fwrite(&buf1, sizeof(unsigned int), 1, fptout);
    }

    printf("Finished Compression!\n");
    printf("Compressed file is called compressed.huf\n");

    fclose(fptout);
    free(in);
}

//decompresses file
void decompress(char *file) {
    printf("Decompressing ...\n");
    FILE *fptdict;
    int freq[256],len_code[256];
    unsigned int code[256];
    unsigned char junk;
    int i,j,total_chars,num_chars=0;
    
    //open dictionary file 
    fptdict = fopen("dictionary", "rb");
    //(format is length space code, with length as int and code as uns char)
    
    //read in file, loop with 3 reads. freq,length,code
    for(i=0;i<256;++i){
        //read to freq array
        fread(&freq[i],sizeof(int),1,fptdict);
        //read space
        fread(&junk,1,1,fptdict);
        //read length
        fread(&len_code[i],sizeof(int),1,fptdict);
        //read space
        fread(&junk,1,1,fptdict);
        //read code
        fread(&code[i],sizeof(unsigned int),1,fptdict);
        //read newline
        fread(&junk,1,1,fptdict);
    }
    //accumulate the count for total chars
    for(i=0;i<256;++i)total_chars += freq[i];
    
    fclose(fptdict);
    
    FILE            *fptcompress, *fptdecompress;
    unsigned int    buf=0,temp=0,buf2=0,temp2=0;
    unsigned char   output;
    int             buflen = 32,k,m,n=1,templen=0,templen2=0,carryover = 0;
    
    //open compressed and new decompressed file
    fptcompress = fopen(file,"rb");
    fptdecompress = fopen("decompressed.huf","wb");
    if(fptcompress==NULL || fptdecompress==NULL)
        printf("Error compress/decompress files\n");
    //start while loop that ends when we reach eof
    while(fread(&buf,sizeof(unsigned int),1,fptcompress)==1){
        for(j=0;j<32;++j){
            //add to templength
            ++templen;
            //shift buffer over to get similar to code read in from dictionary
            temp=buf>>(32-templen);
            //search through "codes" with length as guide
            for(i=0;i<256;++i){
                //if we find a match,
                if(len_code[i] == templen && code[i]==temp) {
                    //output index as unsigned char to decompressed file
                    output = (unsigned char)i;
                    fwrite(&output,1,1,fptdecompress);
                    //if we've written all the chars then quit
                    ++num_chars;
                    if(num_chars==total_chars){
                        printf("Finished Decompression!\n");
                        printf("Decompressed file is called decompressed.huf\n");
                        fclose(fptdecompress);
                        fclose(fptcompress);
                        exit(0);
                    }
                    //reset the buffer and its size
                    buf = buf << len_code[i];
                    buflen -= len_code[i];
                    if(buflen==0) buflen=32;
                    temp=0;templen=0;
                    break;
                }
            }
            //split over two different buffers
            if(j==31 && templen!=0){
                carryover = 0;
                //get the next unsigned int
                if(fread(&buf2,sizeof(unsigned int),1,fptcompress)==1){
                    ++n;
                    for(k=0;k<32;++k){
                        //increment length of temp2
                        ++templen2;
                        //get the new bit into temp so we can compare code
                        temp2 = buf2>>(32-templen2);
                        temp = temp<<1;
                        temp = temp|temp2;
                        ++templen;
                        for(m=0;m<256;++m){
                            if(len_code[m] == templen && code[m]==temp) {
                                //output to file
                                output = (unsigned char)m;
                                fwrite(&output,1,1,fptdecompress);
                                //if we've written all the chars then quit
                                ++num_chars;
                                if(num_chars==total_chars){
                                    printf("Finished Decompression!\n");
                                    printf("Decompressed file is called decompressed.huf\n");
                                    fclose(fptcompress);
                                    fclose(fptdecompress);
                                    exit(0);
                                }
                                //set buf to new buf with the used bits shifted out
                                buf = buf2<<(templen2);
                                //set j to the number of bits lost from buf
                                j=templen2-1;
                                buflen=32-templen2;
                                //clear all of the temps
                                temp=0;temp2=0;templen=0;templen2=0;buf2=0;
                                carryover = 1;
                                break;
                            }
                        }if(carryover==1)break;
                    }
                }
            }                 
        }      
        ++n;                  
    }

    fclose(fptcompress);
    fclose(fptdecompress);
}

//main takes in cmd line args and calls corresponding function
int main(int argc, char *argv[]) { 
    if (argc != 3) {
        printf("Error! Formatting is: \"executable\" [filename] [c/d]\n");
        exit(0);
    }
    if (strcmp(argv[2],"c") != 0) {
        if (strcmp(argv[2],"d") != 0) {
            printf("Error! Formatting is: \"executable\" [filename] [c/d]\n");
            printf("Where c is for compression, d is for decompression\n");
            exit(0);
        }
    }

    if (strcmp(argv[2],"c") == 0) 
        compress(argv[1]);
    if (strcmp(argv[2],"d") == 0)
        decompress(argv[1]);

    return 0; 
}
