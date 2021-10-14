
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Darwin.h"

#define BUFFER_SIZE 1024  /* Maximum length of a line in input file */

/* Uncomment the following line to enable debugging prints 
 * or comment to disable it */
#define DEBUG

#ifdef DEBUG
#define DPRINT(...) fprintf(stderr, __VA_ARGS__);
#else  /* DEBUG */
#define DPRINT(...)
#endif /* DEBUG */
struct ContinentPopulation* deleteFromContinentsTree(struct ContinentPopulation *root,struct Population *node);
struct Species* deleteFromSpecies(struct Species *root,int sid);
void deleteAllPopulations(struct Population *root);
int print_continents();

/**
 * @brief Optional function to initialize data structures that 
 *        need initialization
 *
 * @return 1 on success
 *         0 on failure
 */
int initialize (void)
{
	Species_root=NULL;
    Homo_Sapiens_root=NULL;
    int i=0;
    for(i=0;i<5;i++){
        continents[i]=malloc(sizeof(struct ContinentTree));
        continents[i]->population_root=malloc(sizeof(struct ContinentPopulation));
        continents[i]->sentinel=malloc(sizeof(struct ContinentPopulation));
        continents[i]->sentinel->gid=-1;
        continents[i]->sentinel->lc=NULL;
        continents[i]->sentinel->rc=NULL;
        continents[i]->population_root=continents[i]->sentinel;
    }
	return 1;
}

//inorder gia Population tree
void inorder(struct Population* r){
    
    if(r==NULL)
        return;
    
    inorder(r->lc);
    printf("<%d,%d>",r->gid,r->sid);
    inorder(r->rc);
    
}

void printPopulation(struct Population* r){
   
    if(r==NULL)
        return;

    printPopulation(r->lc);
    printf("<%d>",r->gid);
    printPopulation(r->rc);
    
}
//inorder gia continent tree
void inorderContinents(struct ContinentPopulation* r){
    
    if(r==NULL)
        return;
    inorderContinents(r->lc);
    if(r->gid!=-1)
    printf("<%d>",r->gid);
    inorderContinents(r->rc);
}
void printPop(struct Population*root){
    if(root==NULL)
        return;
    printPop(root->lc);
    printf("<%d,%d>",root->gid,root->cid);
    printPop(root->rc);
}
//postorder gia species
void postorder(struct Species *t){
    if(t==NULL) 
        return;
    postorder(t->lc);
    postorder(t->rc);
    printf("<%d>\n ",t->sid);
    printPop(t->population_root);
    printf("\n");
}
void homoPost(struct HomoSapiens *t){
    if(t==NULL) 
        return;
    homoPost(t->lc);
    homoPost(t->rc);
    
    if(t->lc==NULL&&t->rc==NULL){
      printf("<%d>\n",t->sid);
      printf(" ");
      printPop(t->population_root);
      printf("\n");
    }

}
void MergePrint(struct Species *t){
    if(t==NULL) 
        return;
    MergePrint(t->lc);
    MergePrint(t->rc);
    printf("<%d>\n ",t->sid);
    printPopulation(t->population_root);
    printf("\n");
}
void printSpecies(struct Species *root){
    if(root==NULL)
        return;
    printSpecies(root->lc);
    printSpecies(root->rc);
    printf("<%d>",root->sid);
}
// voithiki gia insert species
struct Species* help_insert(struct Species *r,int sid){
    if(r==NULL){
        struct Species *new=malloc(sizeof(struct Species));
        new->lc=NULL;
        new->rc=NULL;
        new->population_root=NULL;
        new->sid=sid;
        r=new;
        return r;
    }
    if(r->sid==sid){
        return r;
    }
    if(r->sid<sid){
        struct Species * new=malloc(sizeof(struct Species));
        new->lc=r;
        new->population_root=NULL;
        new->rc=NULL;
        new->sid=sid;
        return new;
    }

    if(r->lc!=NULL&&r->lc->sid<sid)
        r->rc=help_insert(r->rc,sid);
    else 
       r->lc=help_insert(r->lc,sid);
  
}









/**
 * @brief insert new species in Species' list with ID <sid>
 *
 * @return 1 on success
 *         0 on failure
 */
int insert_species (int sid)
{
	Species_root=help_insert(Species_root,sid);
    printf("S<%d>\n  ",sid);
    printSpecies(Species_root);
    printf("\nDONE\n");
    return 1;

	return 1;
}


//find the certain species node
struct Species * findnode(struct Species *r,int sid){
    if(r==NULL)
        return NULL;
    if(r->sid==sid)
        return r;
    
    if(r->lc!=NULL&&r->lc->sid<sid)
        findnode(r->rc,sid);
    else
        findnode(r->lc,sid);
   
}
//insert population voithiki
struct Population *ip(struct Population *r,int gid,int sid,int cid,struct Population * parent){
    if(r==NULL){
        struct Population *new=malloc(sizeof(struct Population));
        new->cid=cid;
        new->gid=gid;
        new->sid=sid;
        new->rc=NULL;
        new->lc=NULL;
        new->parent=parent;
        r=new;
        return r;
    }
    if(r->gid==gid)
        return r;
    
    if(r->gid<gid)
        r->rc=ip(r->rc,gid,sid,cid,r);
    else
        r->lc=ip(r->lc,gid,sid,cid,r);

    return r;
}


/**
 * @brief insert new population with ID <gid> in Species' with ID <sid>
 *
 * @return 1 on success
 *         0 on failure
 */
int insert_population(int gid, int sid, int cid){
	struct Species *tmp=findnode(Species_root,sid);
    if(tmp==NULL)
        return 0;
    tmp->population_root=ip(tmp->population_root,gid,sid,cid,NULL);
    printf("G<%d><%d><%d>\n<%d>\n  ",gid,sid,cid,sid);
    printPop(tmp->population_root);
    printf("\nDONE\n");
	return 1;
}



//lowest common ancestor
struct Population *parent(struct Population *r,struct Population *firstChild,struct Population *secondChild){
    if(r->gid<firstChild->gid&&r->gid<secondChild->gid)
        parent(r->rc,firstChild,secondChild);
    else if(r->gid>firstChild->gid&&r->gid>secondChild->gid)
        parent(r->lc,firstChild,secondChild);
    else
        return r;
}
//find certain population node
struct Population *findPopNode(struct Population *r,int gid){
    if(r==NULL)
        return r;
    if(r->gid==gid)
        return r;
    if(r->gid<gid)
        findPopNode(r->rc,gid);
    else
        findPopNode(r->lc,gid);
}


/**
 * @brief find the lowest (earliest) common ancestor of populations with ID <gid1> and <gid2> of species with ID <sid>
 *
 * @return 1 on success
 *         0 on failure
 */
int lowest_common_ancestor(int sid, int gid1, int gid2){
	struct Species *tmp=findnode(Species_root,sid);
    struct Population *first=findPopNode(tmp->population_root,gid1);
    struct Population *second=findPopNode(tmp->population_root,gid2),*commonParent;
    if(first==NULL||second==NULL)
        return 0;
    commonParent=parent(tmp->population_root,first,second);
    printf("L<%d><%d><%d>\n  ",sid,gid1,gid2);
    printf("Lowest Common Ancestor: %d\n",commonParent->gid);
    printf("DONE\n");
	return 1;
}

// count size of population tree
int count_size(struct Population *r){
    if (r==NULL)  
    return 0; 
    else     
    return(count_size(r->lc) + 1 + count_size(r->rc)); 
}

int counter=0;
//copy elements of population tree to given array
void fill_array(struct Population *r,int *array,int *Cid){
    if(r==NULL)
        return;
    fill_array(r->lc,array,Cid);
    array[counter] = r->gid;
    Cid[counter]= r->cid;
    counter++;
    fill_array(r->rc,array,Cid);
}


struct Population * arrayToBST(int *array,int *fcid,int start,int end,int sid){
    if(start>end)
        return NULL;
    int mid=(start+end)/2;
    struct Population *new=malloc(sizeof(struct Population));
    new->gid=array[mid];
    new->cid=fcid[mid];
    new->lc=NULL;
    new->rc=NULL;
    new->sid=sid;
    new->lc=arrayToBST(array,fcid,start,mid-1,sid);
    new->rc=arrayToBST(array,fcid,mid+1,end,sid);

    return new;
}
void merge_arrays(int *a,int*b,int * c,int size1,int size2,int *cid1,int*cid2,int*cid3){
    int i = 0, j = 0, k = 0;
    for (i = 0; i < 11; i++)
    {
        if (j<size1||k<size2)
        {
            if(j==size1){
                c[i]=b[k];
                cid3[i]=cid2[k];
                k++;
            }
            else if(k==size2){
                c[i]=a[j];
                cid3[i]=cid1[j];
                j++;
            }
            else if (a[j] < b[k])
            {
                c[i] = a[j];
                cid3[i]=cid1[j];
                j++;
            }
            else if(a[j]>b[k]){
                c[i]=b[k];
                cid3[i]=cid2[k];
                k++;
            }
        }
    }
}
//helping function to restore all the parent pointers
void restoreParent(struct Population *r,struct Population *prev){
    if(r==NULL)
        return;
    restoreParent(r->lc,r);
    r->parent=prev;
    restoreParent(r->rc,r);
}






/**
 * @brief merge species with IDs <sid1> and <sid2> into new species with ID <sid3>
 *
 * @return 1 on success
 *         0 on failure
 */
int merge_species(int sid1, int sid2, int sid3){
    struct Species *first=findnode(Species_root,sid1),*second=findnode(Species_root,sid2),*third;
    insert_species(sid3);
    third=findnode(Species_root,sid3);

    if(first==NULL||second==NULL)
        return 0;
    int first_size=count_size(first->population_root),second_size=count_size(second->population_root);
    int *firstArr=malloc(first_size*sizeof(int)),*secondArr=malloc(second_size*sizeof(int)),*finalArr=malloc((first_size+second_size)*sizeof(int));
    int *firstCid=malloc(first_size*sizeof(int)),*secondCid=malloc(second_size*sizeof(int)),*finalCid=malloc((first_size+second_size)*sizeof(int));
    counter=0;
    fill_array(first->population_root,firstArr,firstCid);
    counter=0;  
    fill_array(second->population_root,secondArr,secondCid);
    merge_arrays(firstArr,secondArr,finalArr,first_size,second_size,firstCid,secondCid,finalCid);
    third->population_root=arrayToBST(finalArr,finalCid,0,(first_size+second_size)-1,sid3);
    restoreParent(third->population_root,NULL);
    deleteAllPopulations(first->population_root);
    deleteAllPopulations(second->population_root);
    int f=first->sid,s=second->sid;
    Species_root=deleteFromSpecies(Species_root,f);
    Species_root=deleteFromSpecies(Species_root,s);
    printf("M<%d>,<%d>,<%d>\n",sid1,sid2,sid3);
    MergePrint(Species_root);
    return 1;
}

struct ContinentPopulation* addToContinents(struct Population *node,struct ContinentPopulation*rt){//den ginetai na ftasei to O(n) 
    if(rt==continents[node->cid]->sentinel){
        struct ContinentPopulation *new=malloc(sizeof(struct ContinentPopulation));
        new->gid=node->gid;
        new->lc=continents[node->cid]->sentinel;
        new->rc=continents[node->cid]->sentinel;
        return new;
    }
    if(node->gid>rt->gid){
        struct ContinentPopulation *new=malloc(sizeof(struct ContinentPopulation));
        new->gid=node->gid;
        new->lc=rt;
        new->rc=continents[node->cid]->sentinel;
        return new;
    }
    if(node->gid>continents[node->cid]->population_root->gid)
        rt->rc=addToContinents(node,rt->rc);
    else
        rt->lc=addToContinents(node,rt->lc);

    return rt;

}
//preorder population
void preorderPopulation(struct Population *r){//O(n) 
    if(r==NULL)
        return;
    continents[r->cid]->population_root=addToContinents(r,continents[r->cid]->population_root);
    preorderPopulation(r->lc);
    preorderPopulation(r->rc);
}
//preorder for species
void preorderSpecies(struct Species *r){
    if(r==NULL)
        return;

    preorderPopulation(r->population_root);
    preorderSpecies(r->lc);
    preorderSpecies(r->rc);
}




/**
 * @brief Distribute species' population in continents array
 *
 * @return 1 on success
 *         0 on failure
 */
int distribute(){
	preorderSpecies(Species_root);
    int i=0;
    printf("D\n");
    for(i=0;i<5;i++){
        printf("Continent %d :",i);
        inorderContinents(continents[i]->population_root);
        printf("\n");
    }
    printf("DONE\n");
	return 1;
}


//next inorder for Population tree
struct Population* nextInorder(struct Population *root){
    struct Population *tmp=root,*temp=NULL;
    root=root->rc;
    while(root->lc!=NULL)//terma aristera to amesws megalitero tou root
    {
        temp=root;
        root=root->lc;
    }
    tmp->gid=root->gid;
    
    if(root->rc==NULL){
        if(temp!=NULL)
        temp->lc=NULL;
        root=NULL;
        return tmp;
    }

    temp=root->rc;
    root->gid=root->rc->gid;
    root->lc=root->rc->lc;
    root->rc=root->rc->rc;
   // free(temp);
    temp=NULL;
    return tmp;
    
}
//delete from population tree
struct Population* deleteFromPopulation(struct Population *root,int gid){
    if(root==NULL)
        return NULL;
    
    if(root->gid==gid){
        if(root->lc==NULL&&root->rc==NULL){
            continents[root->cid]->population_root=deleteFromContinentsTree(continents[root->cid]->population_root,root);
            
           // free(root);
            root=NULL;
            return NULL;
        }
        else if(root->lc==NULL){
            struct Population *tmp=root->rc;
            continents[root->cid]->population_root=deleteFromContinentsTree(continents[root->cid]->population_root,root);
           
           // free(root);
            root=NULL;
            return tmp;
        }
        else if(root->rc==NULL){
            struct Population*tmp=root->lc;
            continents[root->cid]->population_root=deleteFromContinentsTree(continents[root->cid]->population_root,root); 
           
           // free(root);
            root=NULL;
            return tmp;
        }
        else{
            
            continents[root->cid]->population_root=deleteFromContinentsTree(continents[root->cid]->population_root,root);
    
            root=nextInorder(root);
            return root;
        }
    }

    if(root->gid<gid)
       root->rc= deleteFromPopulation(root->rc,gid);
    else
       root->lc= deleteFromPopulation(root->lc,gid);

    return root;
}
//next inorder for continents tree
struct ContinentPopulation * nextInorderContinent(struct ContinentPopulation *root,int cid){
    struct ContinentPopulation * tmp=root,*temp=NULL;
    root=root->rc;
    while(root->lc->gid!=-1){
        temp=root;
        root=root->lc;
    }   
    tmp->gid=root->gid;
    if(root->rc->gid==-1){
        //free(root);
        tmp->lc=NULL;
        root=NULL;
        return tmp;
    }

    temp=root->rc;
    root->gid=root->rc->gid;
    root->lc=root->rc->lc;
    root->rc=root->rc->rc;
   // free(temp);
    temp=NULL;
    return tmp;
    
}
//delete from continents tree
struct ContinentPopulation* deleteFromContinentsTree(struct ContinentPopulation *root, struct Population *node){
    if(root==continents[node->cid]->sentinel)
        return root;

    if(root->gid==node->gid){
        if(root->lc->gid==-1&&root->rc->gid==-1){
           // free(root);
            root=NULL;
            return continents[node->cid]->sentinel;
        }
        else if(root->lc->gid==-1){
            
            struct ContinentPopulation * tmp=root->rc;
            //free(root);
            root=NULL;
            return tmp;
        }
        else if(root->rc->gid==-1){
            struct ContinentPopulation * tmp=root->lc;
           // free(root);
            root=NULL;
            return tmp;
        }
        else{
            root=nextInorderContinent(root,node->cid);
            return root;
        }
    }



    
    if(root->gid<node->gid)
       root->rc= deleteFromContinentsTree(root->rc,node);
    else
       root->lc= deleteFromContinentsTree(root->lc,node);

    return root; 
}
int getCid(struct Population *root,int gid){
    if(root==NULL)
	return 0;
    if(root->gid==gid)
        return root->cid;
    
    if(root->gid>gid)
        getCid(root->lc,gid);
    else
        getCid(root->rc,gid);
}

/**
 * @brief delete population with ID <gid> from species with ID <sid>
 *
 * @return 1 on success
 *         0 on failure
 */
int delete_population(int gid, int sid){
	struct Species *toDelete=findnode(Species_root,sid);
    
    if(toDelete==NULL)
        return 0;
   
    int cid=getCid(toDelete->population_root,gid);    
 	   
    toDelete->population_root=deleteFromPopulation(toDelete->population_root,gid);
    printf("K<%d><%d>\nSPECIES\n",gid,sid);
    printf("<%d>\n  ",sid);
    printPopulation(toDelete->population_root);
    printf("\nCONTINENTS\n  continent<%d>:",cid);
    inorderContinents(continents[cid]->population_root);
    printf("\nDONE\n");
	
	return 1;
}


//find smallest species node
struct Species *findSmaller(struct Species *root){//O(h) ena monopoati akolouthei
    if(root->lc==NULL&&root->rc==NULL)
        return root;
    
        

    if(root->lc==NULL)
        findSmaller(root->rc);
    else
        findSmaller(root->lc);
}
//delete all the nodes of the given population tree
void deleteAllPopulations(struct Population *root){
    if(root==NULL)
        return; 
    deleteAllPopulations(root->lc);
    deleteAllPopulations(root->rc);
    delete_population(root->gid,root->sid);

}
struct Species* deleteFromSpecies(struct Species *root,int sid){
    if(root==NULL)
        return NULL;
    if(root->sid==sid){
        if(root->rc==NULL&&root->lc==NULL){
            //free(root);
            root=NULL;
            return NULL;
        }
        else if(root->lc==NULL){
            struct Species *tmp=root->rc;
           // free(root);
            root=NULL;
            return tmp;
        }
        else if(root->rc==NULL){
            struct Species *tmp=root->lc;
           // free(root);
            root=NULL;
            return tmp;
        }
        else{
            struct Species *tmp=root->rc,*temp=NULL; 
            while(tmp->lc!=NULL)
                tmp=tmp->lc;
            
            tmp->lc=root->lc;
            root->lc=tmp->lc;
            
            root->sid=root->rc->sid;
            temp=root->rc;
            root->rc=root->rc->rc;
           // free(temp);
            temp=NULL;
            
            return root;
        }

    }

    
    if(root->lc!=NULL&&root->lc->sid<sid)
        root->rc=deleteFromSpecies(root->rc,sid);
    else if(root->lc!=NULL)
        root->lc=deleteFromSpecies(root->lc,sid);
    else
        root->rc=deleteFromSpecies(root->rc,sid);
    
    return root;
}

/**
 * @brief delete the species with lowest id and its populations
 *
 * @return 1 on success
 *         0 on failure
 */
int delete_species(){
	struct Species *smallest=findSmaller(Species_root);
    if(smallest->population_root!=NULL);
        deleteAllPopulations(smallest->population_root);
    
    Species_root=deleteFromSpecies(Species_root,smallest->sid);
    
    printf("F\nSPECIES\n");
    MergePrint(Species_root);
    printf("CONTINENTS\n");
    print_continents();
	return 1;
}

int Max(int a,int b){
    if(a>b)
        return a;
    return b;
}
int Min(int a,int b){
    if(a<b)
        return a;
    return b;
}

struct HomoSapiens* insertToHomoSapiens(struct HomoSapiens * root, struct Species *species){
   
    if(root==NULL){
        struct HomoSapiens *new=malloc(sizeof(struct HomoSapiens));
        new->lc=NULL;
        new->rc=NULL;
        new->next=NULL;
        new->sid=species->sid;
        new->population_root=species->population_root;
        root=new;
        return root;
    }
    
    if(root->rc==NULL&&root->lc==NULL){
        struct HomoSapiens *new ,*new1;
        new= malloc(sizeof(struct HomoSapiens));
        new1=malloc(sizeof(struct HomoSapiens));
        root->lc=new;
        root->rc=new1;
        new->lc=NULL;
        new->rc=NULL;
        new->next=NULL;
        new1->rc=NULL;
        new1->lc=NULL;
        new1->next=NULL;
        new->sid=Min(root->sid,species->sid);
        new1->sid=Max(root->sid,species->sid);
        if(root->sid==Min(root->sid,species->sid)){
            new->population_root=root->population_root;
            new1->population_root=species->population_root;
        }
        else{
            new->population_root=species->population_root;
            new1->population_root=root->population_root;
            root->population_root=new->population_root;
        }
        new1->next=root->next;
        root->next=new1;
        new->next=new1;
        root->sid=new->sid;
        return root;
    }
    if(root->sid>=species->sid)
         insertToHomoSapiens(root->lc,species);
    else
         insertToHomoSapiens(root->rc,species);
    return root;
}
void create_homo_tree(struct Species *root){
    if(root==NULL)
        return;
    Homo_Sapiens_root=insertToHomoSapiens(Homo_Sapiens_root,root);    
    create_homo_tree(root->lc);
    create_homo_tree(root->rc);
}

void restoreNext(struct HomoSapiens *root){
    while(root->lc!=NULL)
        root=root->lc;
    struct HomoSapiens *tmp=root;
    while(tmp->next!=NULL){
        tmp=tmp->next;
        while(tmp->lc!=NULL)
            tmp=tmp->lc;
        root->next=tmp;
        root=root->next;
    }
}


void inorder_homo(struct HomoSapiens *root){
    if(root==NULL)
        return;
    inorder_homo(root->lc);
    if(root->rc==NULL&&root->lc==NULL){
        if(root->population_root!=NULL)
        printf("[");
        inorder(root->population_root);
        if(root->population_root!=NULL)
        printf("]");
    }
    
    inorder_homo(root->rc);
    
}



/**
 * @brief Remaining species evolve into homo sapiens.
 *
 * @return 1 on success
 *         0 on failure
 */
int evolution(){

	if(Species_root==NULL)
        return 0;
    printf("E\n  ");
    printf("Homo Sapiens:");
    create_homo_tree(Species_root);
    inorder_homo(Homo_Sapiens_root);
    printf("\nDONE\n");
    restoreNext(Homo_Sapiens_root);
    while(Species_root!=NULL)
        Species_root=deleteFromSpecies(Species_root,Species_root->sid);
	return 1;
}

/**
 * @brief Gather statistics on species from homo_sapiens tree
 *
 * @return 1 on success
 *         0 on failure
 */
int species_statistics(){
	struct HomoSapiens *root=Homo_Sapiens_root;
    if(root==NULL)
        return 0;
    printf("N\n  Homo Sapiens:");
    while(root->lc!=NULL)
        root=root->lc;
    int sum=0;
    
    while(root!=NULL){
        if(root->next!=NULL)
            printf("<%d>,",root->sid);
        else
             printf("<%d>\n  ",root->sid);
        sum++;
        root=root->next;
    }
    printf("Homo Sapiens species:<%d>\n",sum);
    printf("DONE\n");
	return 1;
}

struct HomoSapiens *findHomoNode(struct HomoSapiens *root,int sid){
    if(root==NULL)
        return NULL;
    if(root->sid==sid&&root->lc==NULL&&root->rc==NULL)
        return root;
    
    if(root->sid>=sid)
        findHomoNode(root->lc,sid);
    else
        findHomoNode(root->rc,sid);
}

/**
 * @brief Gather statistics on population belonging to species with ID <sid> from homo_sapiens tree
 *
 * @return 1 on success
 *         0 on failure
 */
int population_statistics(int sid){
	struct HomoSapiens *root=findHomoNode(Homo_Sapiens_root,sid);
    if(root==NULL)
        return 0;
    printf("J<%d>\n  ",sid);
    int size=count_size(root->population_root);
    printf("Homo Sapiens species:<%d>\nDONE",size);
	return 1;
}

/**
 * @brief Print species' leaves list
 *
 * @return 1 on success
 *         0 on failure
 */
int print_species(){
	printf("P\n  ");
    if(Species_root==NULL){ // an ginei evolution alliws tupwnei kanonika to species tree
        struct HomoSapiens *root=Homo_Sapiens_root;
         while(root->lc!=NULL)
            root=root->lc;
    
    
        while(root!=NULL){
            if(root->next!=NULL)
            printf("<%d>,",root->sid);
            else
            printf("<%d>",root->sid);
            root=root->next;
        }
        
    }  
    
    printSpecies(Species_root);
    printf("\nDONE\n");
	return 1;
}

/**
 * @brief Print species' tree with population trees
 *
 * @return 1 on success
 *         0 on failure
 */
int print_populations(){
	printf("X\n");
    homoPost(Homo_Sapiens_root);
    printf("DONE\n");
    return 1;
}

/**
 * @brief Print continents array with each continent's population tree
 *
 * @return 1 on success
 *         0 on failure
 */
int print_continents(){
	int i=0;
    for(i=0;i<5;i++){
        printf("Continent %d :",i);
        inorderContinents(continents[i]->population_root);
        printf("\n");
    }
    printf("DONE\n");
    return 1;
}

/**
 * @brief Print homo_sapiens tree
 *
 * @return 1 on success
 *         0 on failure
 */
int print_homo_sapiens(){
	printf("W\n Homo Sapiens:");
    inorder_homo(Homo_Sapiens_root);
    printf("\nDONE\n");
    return 1;
}

/**
 * @brief Free resources
 *
 * @return 1 on success
 *         0 on failure
 */
int free_all(void)
{
	return 1;
}


/**
 * @brief The main function
 *
 * @param argc Number of arguments
 * @param argv Argument vector
 *
 * @return 0 on success
 *         1 on failure
 */
int main(int argc, char** argv)
{
	FILE *fin = NULL;
	char buff[BUFFER_SIZE], event;

	/* Check command buff arguments */
	if (argc != 2) {
		fprintf(stderr, "Usage: %s <input_file> \n", argv[0]);
		return EXIT_FAILURE;
	}

	/* Open input file */
	if (( fin = fopen(argv[1], "r") ) == NULL ) {
		fprintf(stderr, "\n Could not open file: %s\n", argv[1]);
		perror("Opening test file\n");
		return EXIT_FAILURE;
	}

	/* Initializations */
	initialize();

	/* Read input file buff-by-buff and handle the events */
	while (fgets(buff, BUFFER_SIZE, fin)) {

		DPRINT("\n>>> Event: %s", buff);

		switch(buff[0]) {

			/* Comment */
			case '#':
				break;

				/* Insert Species
				 * S <sid> */
			case 'S':
				{
					int sid;

					sscanf(buff, "%c %d", &event, &sid);
					DPRINT("%c %d\n", event, sid);

					if (insert_species(sid)) {
						DPRINT("%c %d succeeded\n", event, sid);
					} else {
						fprintf(stderr, "%c %d failed\n", event, sid);
					}

					break;
				}

				/* Insert population
				 * G <gid> <sid> <cid> */
			case 'G':
				{
					int gid, sid, cid;

					sscanf(buff, "%c %d %d %d", &event, &gid, &sid, &cid);
					DPRINT("%c %d %d %d\n", event, gid, sid, cid);

					if (insert_population(gid, sid, cid)) {
						DPRINT("%c %d %d %d succeeded\n", event, gid, sid, cid);
					} else {
						fprintf(stderr, "%c %d %d %d failed\n", event, gid, sid, cid);
					}

					break;
				}

				/* Lowest Common Ancestor
				 * L <sid> <gid1> <gid2> */
			 case 'L':
				{
					int sid, gid1, gid2;

					sscanf(buff, "%c %d %d %d", &event, &sid, &gid1, &gid2);
					DPRINT("%c %d %d %d\n", event, sid, gid1, gid2);

					if (lowest_common_ancestor(sid, gid1, gid2)) {
						DPRINT("%c %d %d %d succeeded\n", event, sid, gid1, gid2);
					} else {
						fprintf(stderr, "%c %d %d %d failed\n", event, sid, gid1, gid2);
					}

					break;
				}

				/* Merge species
				 * M <sid1> <sid2> <sid3> */
			case 'M':
				{
					int sid1, sid2, sid3;

					sscanf(buff, "%c %d %d %d", &event, &sid1, &sid2, &sid3);
					DPRINT("%c %d %d %d\n", event, sid1, sid2, sid3);

					if (merge_species(sid1, sid2, sid3)) {
						DPRINT("%c %d %d %d succeeded\n", event, sid1, sid2, sid3);
					} else {
						fprintf(stderr, "%c %d %d %d failed\n", event, sid1, sid2, sid3);
					}

					break;
				}

				/* Distribute species
				 * D */
			case 'D':
				{
					sscanf(buff, "%c", &event);
					DPRINT("%c\n", event);

					if (distribute()) {
						DPRINT("%c succeeded\n", event);
					} else {
						fprintf(stderr, "%c failed\n", event);
					}

					break;
				}

				/* Delete population
				 * K <gid> <sid> */
			case 'K':
				{
					int gid, sid;

					sscanf(buff, "%c %d %d", &event, &gid, &sid);
					DPRINT("%c %d %d\n", event, gid, sid);

					if (delete_population(gid, sid)) {
						DPRINT("%c %d %d succeeded\n", event, gid, sid);
					} else {
						fprintf(stderr, "%c %d %d failed\n", event, gid, sid);
					}

					break;
				}

				/* Delete species
				 * F */
			case 'F':
				{
					sscanf(buff, "%c", &event);
					DPRINT("%c \n", event);

					if (delete_species()) {
						DPRINT("%c succeeded\n", event);
					} else {
						fprintf(stderr, "%c failed\n", event);
					}

					break;
				}

				/* Evolution to homo sapiens
				 * E */
			case 'E':
				{
					sscanf(buff, "%c", &event);
					DPRINT("%c\n", event);

					if (evolution()) {
						DPRINT("%c succeeded\n", event);
					} else {
						fprintf(stderr, "%c failed\n", event);
					}

					break;
				}

				/* Species' statistics
				 * N */
			case 'N':
				{
					sscanf(buff, "%c", &event);
					DPRINT("%c\n", event);

					if (species_statistics()) {
						DPRINT("%c succeeded\n", event);
					} else {
						fprintf(stderr, "%c failed\n", event);
					}

					break;
				}

				/* Population statistics
				 * J <sid> */
			case 'J':
				{
					int sid;

					sscanf(buff, "%c %d", &event, &sid);
					DPRINT("%c %d\n", event, sid);

					if (population_statistics(sid)) {
						DPRINT("%c %d succeeded\n", event, sid);
					} else {
						fprintf(stderr, "%c %d failed\n", event, sid);
					}

					break;
				}

				/* Print species
				 * P */
			case 'P':
				{
					sscanf(buff, "%c", &event);
					DPRINT("%c\n", event);

					if (print_species()) {
						DPRINT("%c succeeded\n", event);
					} else {
						fprintf(stderr, "%c failed\n", event);
					}

					break;
				}

			/* Print populations
				 * X */
			case 'X':
				{
					sscanf(buff, "%c", &event);
					DPRINT("%c\n", event);

					if (print_populations()) {
						DPRINT("%c succeeded\n", event);
					} else {
						fprintf(stderr, "%c failed\n", event);
					}

					break;
				}

			/* Print continents
				 * C */
			case 'C':
				{
					sscanf(buff, "%c", &event);
					DPRINT("%c\n", event);

					if (print_continents()) {
						DPRINT("%c succeeded\n", event);
					} else {
						fprintf(stderr, "%c failed\n", event);
					}

					break;
				}

			/* Print homo sapiens
				 * W */
			case 'W':
				{
					sscanf(buff, "%c", &event);
					DPRINT("%c\n", event);

					if (print_homo_sapiens()) {
						DPRINT("%c succeeded\n", event);
					} else {
						fprintf(stderr, "%c failed\n", event);
					}

					break;
				}

				/* Empty line */
			case '\n':
				break;

				/* Ignore everything else */
			default:
				DPRINT("Ignoring buff: %s \n", buff);

				break;
		}
	}

	free_all();
	return (EXIT_SUCCESS);
}
