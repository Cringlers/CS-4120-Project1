#include<iostream>
#include<fstream>
#include<iomanip>
#include<chrono>
#include "number_data.h"

using namespace std;
using namespace chrono;

void BubbleSort(int list[], int Size, string listType) {
    int numComps = 0;
    int temp = 0;
    bool Swap = false;
    const auto t1 = high_resolution_clock::now();

    for (int i = 0; i < Size - 1; i++) {
        for (int j = 0; j < Size - i - 1; j++) {
            if (list[j] > list[j + 1]) {
                temp = list[j + 1];
                list[j + 1] = list[j];
                list[j] = temp;

                Swap = true;
            }
            numComps++;
        }

        if (Swap != true) {
            break;
        }
    }

    const auto t2 = high_resolution_clock::now();
    const auto d = duration_cast<duration<double>> (t2 - t1).count();

    cout << "# of comparisons for "  << Size << " item list in " << listType << " order: " << numComps << endl;
    cout << "Time taken: " << d * 1000000 << " microseconds" << endl;
    cout << endl;

    return;
}

void InsertionSort(int list[], int Size, string listType) {
    int numComps = 0;
    const auto t1 = high_resolution_clock::now();
    
    for (int i = 1; i < Size; i++) {
        int curr = list[i];
        int j = i - 1;

        numComps++;
        while (j >= 0 && list[j] > curr) {
            list[j + 1] = list[j];
            j--;

            numComps++;
        }
    }

    const auto t2 = high_resolution_clock::now();
    const auto d = duration_cast<duration<double>> (t2 - t1).count();
    
    cout << "# of comparisons for " << Size << " item list in " << listType << " order: " << numComps << endl;
    cout << "Time taken: " << d * 1000000 << " microseconds" << endl;
    cout << endl;

    return;
}

int QSRcomps = 0;

int PartitionR(int list[], int low, int high) {
    int temp = 0;
    int pivot = list[rand() % (high + 1)];
    int i = low - 1;

    for (int j = low; j <= high; j++) {
        QSRcomps++;

        if (list[j] < pivot) {
            i++;
            temp = list[j];
            list[j] = list[i];
            list[i] = temp;
        }
    }

    temp = list[i + 1];
    list[i + 1] = list[high];
    list[high] = temp;
    return(i + 1);
}

void QuickSortR(int list[], int low, int high) {
    if (low < high) {

        int parIndex = PartitionR(list, low, high);

        QuickSortR(list, low, parIndex - 1);
        QuickSortR(list, parIndex + 1, high);
    }

    return;
}

int QSMcomps = 0;

int PartitionM(int list[], int low, int high) {
    int temp = 0;
    int first = list[low];
    int last = list[high];
    int mid = list[low + ((high - low) / 2)];
    int i = low - 1;
    int pivot;

    if ((first > last) ^ (first > mid)) {
        pivot = first;
    }
    else if ((last > first) ^ (last > mid)) {
        pivot = last;
    }
    else {
        pivot = mid;
    }

    for (int j = low; j <= high; j++) {
        QSMcomps++;

        if (list[j] < pivot) {
            i++;
            temp = list[j];
            list[j] = list[i];
            list[i] = temp;
        }
    }

    temp = list[i + 1];
    list[i + 1] = list[high];
    list[high] = temp;
    return(i + 1);
}

void QuickSortM(int list[], int low, int high) {
    if (low < high) {

        int parIndex = PartitionM(list, low, high);

        QuickSortM(list, low, parIndex - 1);
        QuickSortM(list, parIndex + 1, high);
    }

    return;
}

void QuickSort1(int list[], int Size, string listType) {
    const auto t1 = high_resolution_clock::now();

    QuickSortR(list, 0, Size - 1);

    const auto t2 = high_resolution_clock::now();
    const auto d = duration_cast<duration<double>> (t2 - t1).count();

    cout << "# of comparisons for " << Size << " item list in " << listType << " order: " << QSRcomps << endl;
    cout << "Time taken: " << d * 1000000 << " microseconds" << endl;
    cout << endl;

    QSRcomps = 0;

    return;
}

void QuickSort2(int list[], int Size, string listType) {
    const auto t1 = high_resolution_clock::now();

    QuickSortM(list, 0, Size - 1);

    const auto t2 = high_resolution_clock::now();
    const auto d = duration_cast<duration<double>> (t2 - t1).count();

    cout << "# of comparisons for " << Size << " item list in " << listType << " order: " << QSMcomps << endl;
    cout << "Time taken: " << d * 1000000 << " microseconds" << endl;
    cout << endl;

    QSMcomps = 0;

    return;
}

int heapComps = 0;

void Heapify(int list[], int Size, int i) {
    int max = i;
    int left = (i * 2) + 1;
    int right = (i * 2) + 2;

    heapComps++;
    if ((left < Size) && (list[left] > list[max])) {
        max = left;
    }

    heapComps++;
    if ((right < Size) && (list[right] > list[max])) {
        max = right;
    }

    if (max != i) {
        int temp = list[max];
        list[max] = list[i];
        list[i] = temp;

        Heapify(list, Size, max);
    }

    return;
}

void HeapSort(int list[], int Size, string listType) {
    const auto t1 = high_resolution_clock::now();
    
    for (int i = (Size / 2) - 1; i >= 0; i--) {
        Heapify(list, Size, i);
    }

    for (int i = Size - 1; i > 0; i--) {
        int temp = list[0];
        list[0] = list[i];
        list[i] = temp;

        Heapify(list, i, 0);
    }

    const auto t2 = high_resolution_clock::now();
    const auto d = duration_cast<duration<double>> (t2 - t1).count();

    cout << "# of comparisons for " << Size << " item list in " << listType << " order: " << heapComps << endl;
    cout << "Time taken: " << d * 1000000 << " microseconds" << endl;
    cout << endl;

    heapComps = 0;

    return;
}

int mergeComps = 0;

void Merge(int Tlist[], int low, int mid1, int mid2, int high, int list[]) {
    int i = low;
    int j = mid1;
    int k = mid2;
    int l = low;

    while ((i < mid1) && (j < mid2) && (k < high)) {
        mergeComps++;

        if (Tlist[i] < Tlist[j]) {
            mergeComps++;

            if (Tlist[i] < Tlist[k]) {
                list[l++] = Tlist[i++];
            }
            else {
                list[l++] = Tlist[k++];
            }
        }
        else {
            mergeComps++;

            if (Tlist[j] < Tlist[k]) {
                list[l++] = Tlist[j++];
            }
            else {
                list[l++] = Tlist[k++];
            }
        }
    }

    while ((i < mid1) && (j < mid2)) {
        mergeComps++;

        if (Tlist[i] < Tlist[j]) {
            list[l++] = Tlist[i++];
        }
        else {
            list[l++] = Tlist[j++];
        }
    }

    while ((j < mid2) && (k < high)) {
        mergeComps++;

        if (Tlist[j] < Tlist[k]) {
            list[l++] = Tlist[j++];
        }
        else {
            list[l++] = Tlist[k++];
        }
    }

    while ((i < mid1) && (k < high)) {
        mergeComps++;

        if (Tlist[i] < Tlist[k]) {
            list[l++] = Tlist[i++];
        }
        else {
            list[l++] = Tlist[k++];
        }
    }

    while (i < mid1) {
        list[l++] = Tlist[i++];
    }

    while (j < mid2) {
        list[l++] = Tlist[j++];
    }

    while (k < high) {
        list[l++] = Tlist[k++];
    }

    return;
}

void MergeSort3Way(int Tlist[], int low, int high, int list[]) {
    if (high - low < 2) {
        return;
    }

    int mid1 = low + ((high - low) / 3);
    int mid2 = low + 2 * ((high - low) / 3) + 1;

    MergeSort3Way(list, low, mid1, Tlist);
    MergeSort3Way(list, mid1, mid2, Tlist);
    MergeSort3Way(list, mid2, high, Tlist);

    Merge(list, low, mid1, mid2, high, Tlist);

    return;
}

void MergeSort(int list[], int Size, string listType) {
    const auto t1 = high_resolution_clock::now();
    
    int* Tlist = new int[Size];

    for (int i = 0; i < Size; i++) {
        Tlist[i] = list[i];
    }

    MergeSort3Way(Tlist, 0, Size, list);

    for (int i = 0; i < Size; i++) {
        list[i] = Tlist[i];
    }

    delete[] Tlist;

    const auto t2 = high_resolution_clock::now();
    const auto d = duration_cast<duration<double>> (t2 - t1).count();

    cout << "# of comparisons for " << Size << " item list in " << listType << " order: " << mergeComps << endl;
    cout << "Time taken: " << d * 1000000 << " microseconds" << endl;
    cout << endl;

    mergeComps = 0;

    return;
}

const int SIZE = 8000;

int main() {
    //----------Create lists----------//
    int* randList = new int[SIZE];
    int* ascList = new int[SIZE];
    int* descList = new int[SIZE];

    int* rand2 = new int[SIZE / 2];
    int* rand3 = new int[SIZE / 4];
    int* rand4 = new int[SIZE / 8];

    int* asc2 = new int[SIZE / 2];
    int* asc3 = new int[SIZE / 4];
    int* asc4 = new int[SIZE / 8];

    int* desc2 = new int[SIZE / 2];
    int* desc3 = new int[SIZE / 4];
    int* desc4 = new int[SIZE / 8];

    fill_random_numbers(randList, SIZE);
    fill_random_numbers(rand2, SIZE / 2);
    fill_random_numbers(rand3, SIZE / 4);
    fill_random_numbers(rand4, SIZE / 8);

    fill_sorted_numbers(ascList, descList, SIZE);
    fill_sorted_numbers(asc2, desc2, SIZE / 2);
    fill_sorted_numbers(asc3, desc3, SIZE / 4);
    fill_sorted_numbers(asc4, desc4, SIZE / 8);

    ofstream dataFile;
    dataFile.open("dataRand.txt");
    for (int i = 0; i < SIZE; i++) {
        dataFile << randList[i] << endl;
    }
    dataFile.close();

    dataFile.open("dataAsc.txt");
    for (int i = 0; i < SIZE; i++) {
        dataFile << ascList[i] << endl;
    }
    dataFile.close();

    dataFile.open("dataDesc.txt");
    for (int i = 0; i < SIZE; i++) {
        dataFile << descList[i] << endl;
    }
    dataFile.close();

    //----------Bubble Sort Calls----------//
    cout << "Bubble Sort:" << endl;
    cout << setw(60) << setfill('-') << "" << endl;

    BubbleSort(randList, SIZE, "random");
    BubbleSort(rand2, SIZE / 2, "random");
    BubbleSort(rand3, SIZE / 4, "random");
    BubbleSort(rand4, SIZE / 8, "random");

    BubbleSort(ascList, SIZE, "ascending");
    BubbleSort(asc2, SIZE / 2, "ascending");
    BubbleSort(asc3, SIZE / 4, "ascending");
    BubbleSort(asc4, SIZE / 8, "ascending");

    BubbleSort(descList, SIZE, "descending");
    BubbleSort(desc2, SIZE / 2, "descending");
    BubbleSort(desc3, SIZE / 4, "descending");
    BubbleSort(desc4, SIZE / 8, "descending");

    cout << endl;

    //----------Refill lists with original data----------//
    ifstream file;
    //-----Random Lists-----//
    int num;
    file.open("dataRand.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        randList[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        rand2[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        rand3[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        rand4[i] = num;
    }
    file.close();

    //-----Ascending Lists-----//
    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        ascList[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        asc2[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        asc3[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        asc4[i] = num;
    }
    file.close();

    //-----Descending Lists-----//
    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        descList[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        desc2[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        desc3[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        desc4[i] = num;
    }
    file.close();

    //----------Insertion Sort Calls----------//
    cout << "Insertion Sort:" << endl;
    cout << setw(60) << setfill('-') << "" << endl;

    InsertionSort(randList, SIZE, "random");
    InsertionSort(rand2, SIZE / 2, "random");
    InsertionSort(rand3, SIZE / 4, "random");
    InsertionSort(rand4, SIZE / 8, "random");

    InsertionSort(ascList, SIZE, "ascending");
    InsertionSort(asc2, SIZE / 2, "ascending");
    InsertionSort(asc3, SIZE / 4, "ascending");
    InsertionSort(asc4, SIZE / 8, "ascending");

    InsertionSort(descList, SIZE, "descending");
    InsertionSort(desc2, SIZE / 2, "descending");
    InsertionSort(desc3, SIZE / 4, "descending");
    InsertionSort(desc4, SIZE / 8, "descending");

    cout << endl;

    //----------Refill lists with original data----------//
    //-----Random Lists-----//
    file.open("dataRand.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        randList[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        rand2[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        rand3[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        rand4[i] = num;
    }
    file.close();

    //-----Ascending Lists-----//
    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        ascList[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        asc2[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        asc3[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        asc4[i] = num;
    }
    file.close();

    //-----Descending Lists-----//
    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        descList[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        desc2[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        desc3[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        desc4[i] = num;
    }
    file.close();

    //----------Quick Sort R Calls----------//
    cout << "Quick Sort R:" << endl;
    cout << setw(60) << setfill('-') << "" << endl;

    QuickSort1(randList, SIZE, "random");
    QuickSort1(rand2, SIZE / 2, "random");
    QuickSort1(rand3, SIZE / 4, "random");
    QuickSort1(rand4, SIZE / 8, "random");

    QuickSort1(ascList, SIZE, "ascending");
    QuickSort1(asc2, SIZE / 2, "ascending");
    QuickSort1(asc3, SIZE / 4, "ascending");
    QuickSort1(asc4, SIZE / 8, "ascending");

    QuickSort1(descList, SIZE, "descending");
    QuickSort1(desc2, SIZE / 2, "descending");
    QuickSort1(desc3, SIZE / 4, "descending");
    QuickSort1(desc4, SIZE / 8, "descending");

    cout << endl;

    //----------Refill lists with original data----------//
    //-----Random Lists-----//
    file.open("dataRand.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        randList[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        rand2[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        rand3[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        rand4[i] = num;
    }
    file.close();

    //-----Ascending Lists-----//
    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        ascList[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        asc2[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        asc3[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        asc4[i] = num;
    }
    file.close();

    //-----Descending Lists-----//
    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        descList[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        desc2[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        desc3[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        desc4[i] = num;
    }
    file.close();

    //----------Quick Sort M Calls----------//
    cout << "Quick Sort M:" << endl;
    cout << setw(60) << setfill('-') << "" << endl;

    QuickSort2(randList, SIZE, "random");
    QuickSort2(rand2, SIZE / 2, "random");
    QuickSort2(rand3, SIZE / 4, "random");
    QuickSort2(rand4, SIZE / 8, "random");

    QuickSort2(ascList, SIZE, "ascending");
    QuickSort2(asc2, SIZE / 2, "ascending");
    QuickSort2(asc3, SIZE / 4, "ascending");
    QuickSort2(asc4, SIZE / 8, "ascending");

    QuickSort2(descList, SIZE, "descending");
    QuickSort2(desc2, SIZE / 2, "descending");
    QuickSort2(desc3, SIZE / 4, "descending");
    QuickSort2(desc4, SIZE / 8, "descending");

    cout << endl;

    //----------Refill lists with original data----------//
    //-----Random Lists-----//
    file.open("dataRand.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        randList[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        rand2[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        rand3[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        rand4[i] = num;
    }
    file.close();

    //-----Ascending Lists-----//
    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        ascList[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        asc2[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        asc3[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        asc4[i] = num;
    }
    file.close();

    //-----Descending Lists-----//
    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        descList[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        desc2[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        desc3[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        desc4[i] = num;
    }
    file.close();

    //----------Heap Sort Calls----------//
    cout << "Heap Sort:" << endl;
    cout << setw(60) << setfill('-') << "" << endl;

    HeapSort(randList, SIZE, "random");
    HeapSort(rand2, SIZE / 2, "random");
    HeapSort(rand3, SIZE / 4, "random");
    HeapSort(rand4, SIZE / 8, "random");

    HeapSort(ascList, SIZE, "ascending");
    HeapSort(asc2, SIZE / 2, "ascending");
    HeapSort(asc3, SIZE / 4, "ascending");
    HeapSort(asc4, SIZE / 8, "ascending");

    HeapSort(descList, SIZE, "descending");
    HeapSort(desc2, SIZE / 2, "descending");
    HeapSort(desc3, SIZE / 4, "descending");
    HeapSort(desc4, SIZE / 8, "descending");

    cout << endl;

    //----------Refill lists with original data----------//
    //-----Random Lists-----//
    file.open("dataRand.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        randList[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        rand2[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        rand3[i] = num;
    }
    file.close();

    file.open("dataRand.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        rand4[i] = num;
    }
    file.close();

    //-----Ascending Lists-----//
    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        ascList[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        asc2[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        asc3[i] = num;
    }
    file.close();

    file.open("dataAsc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        asc4[i] = num;
    }
    file.close();

    //-----Descending Lists-----//
    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE; i++) {
        file >> num;
        descList[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 2; i++) {
        file >> num;
        desc2[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 4; i++) {
        file >> num;
        desc3[i] = num;
    }
    file.close();

    file.open("dataDesc.txt");
    for (int i = 0; i < SIZE / 8; i++) {
        file >> num;
        desc4[i] = num;
    }
    file.close();

    //----------3 Way Merge Sort Calls----------//
    cout << "3 Way Merge Sort:" << endl;
    cout << setw(60) << setfill('-') << "" << endl;

    MergeSort(randList, SIZE, "random");
    MergeSort(rand2, SIZE / 2, "random");
    MergeSort(rand3, SIZE / 4, "random");
    MergeSort(rand4, SIZE / 8, "random");

    MergeSort(ascList, SIZE, "ascending");
    MergeSort(asc2, SIZE / 2, "ascending");
    MergeSort(asc3, SIZE / 4, "ascending");
    MergeSort(asc4, SIZE / 8, "ascending");

    MergeSort(descList, SIZE, "descending");
    MergeSort(desc2, SIZE / 2, "descending");
    MergeSort(desc3, SIZE / 4, "descending");
    MergeSort(desc4, SIZE / 8, "descending");

    return 0;
}