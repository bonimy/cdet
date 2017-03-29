//Calculate median after sorting with HEAPSORT algorithm.
//Also calculate the 15.87% and 84.13% quantiles(q1 and q2).

#include "__base.h"
#include "medsort.h"
#include "sort.h"

int partition(float* list, int left, int right, int pivot);
int select(float* list, int left, int right, int k);

float medsort(float* a, int n, float* q1, float* q2)
{
    if (n <= 1)
        return *a;
    float* as = (float*)malloc(n * sizeof(float));
    memcpy(as, a, n * sizeof(float));

    float amed = medsortdirect(as, n, q1, q2);

    free(as);
    return amed;
}

float medsortdirect(float* a, int n, float* q1, float* q2)
{
    if (n <= 1)
        return *a;

    float amed;
    if (n % 2 == 1)
    {
        amed = a[select(a, 0, n - 1, (n - 1) / 2)];
    }
    else
    {
        amed = a[select(a, 0, n - 1, n / 2 - 1)];
        amed += a[select(a, 0, n - 1, n / 2)];
        amed /= 2;
    }

    if (q1 != 0)
        *q1 = a[select(a, 0, n - 1, max(lroundf(0.1587f*n), 1) - 1)];

    if (q2 != 0)
        *q2 = a[select(a, 0, n - 1, lroundf(0.8413f*n) - 1)];

    return amed;
}

int partition(float* list, int left, int right, int pivot)
{
    float value = list[pivot];
    list[pivot] = list[right];
    list[right] = value;

    int index = left;
    for (int i = left; i < right; i++)
    {
        if (list[i] < value)
        {
            float temp = list[index];
            list[index] = list[i];
            list[i] = temp;
            index++;
        }
    }
    value = list[index];
    list[index] = list[right];
    list[right] = value;
    return index;
}

int select(float* list, int left, int right, int k)
{
    while (1)
    {
        if (left == right)
            return left;

        int med = left + (right - left) / 2;
        int max = list[left] > list[right] ? left : right;
        max = list[max] > list[med] ? max : med;
        int min = list[left] < list[right] ? left : right;
        min = list[min] < list[med] ? min : med;

        int pivot = med;
        if (list[left] > list[min] && list[left] < list[max])
            pivot = left;
        else if (list[right] > list[min] && list[right] < list[max])
            pivot = right;

        pivot = partition(list, left, right, pivot);

        if (k == pivot)
            return k;
        else if (k < pivot)
            right = pivot - 1;
        else
            left = pivot + 1;
    }
}