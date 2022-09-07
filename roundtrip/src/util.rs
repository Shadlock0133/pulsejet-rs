use std::mem::MaybeUninit;

#[warn(unsafe_op_in_unsafe_fn)]
pub const unsafe fn slice_assume_init_ref<T>(slice: &[MaybeUninit<T>]) -> &[T] {
    // SAFETY: casting `slice` to a `*const [T]` is safe since the caller guarantees that
    // `slice` is initialized, and `MaybeUninit` is guaranteed to have the same layout as `T`.
    // The pointer obtained is valid since it refers to memory owned by `slice` which is a
    // reference and thus guaranteed to be valid for reads.
    unsafe { &*(slice as *const [MaybeUninit<T>] as *const [T]) }
}

pub struct ArrayChunks<I: Iterator, const N: usize> {
    iterator: I,
    remainder: Option<([MaybeUninit<I::Item>; N], usize)>,
}

impl<I: Iterator, const N: usize> ArrayChunks<I, N> {
    pub fn remainder(&self) -> Option<&[I::Item]> {
        self.remainder
            .as_ref()
            .map(|(a, i)| &a[..*i])
            .map(|slice| unsafe { slice_assume_init_ref(slice) })
    }
}

impl<I: Iterator, const N: usize> Iterator for ArrayChunks<I, N> {
    type Item = [I::Item; N];

    fn next(&mut self) -> Option<Self::Item> {
        let mut array: [MaybeUninit<I::Item>; N] =
            unsafe { MaybeUninit::uninit().assume_init() };
        for i in 0..N {
            match self.iterator.next() {
                Some(item) => {
                    array[i].write(item);
                }
                None => {
                    self.remainder = Some((array, i));
                    return None;
                }
            }
        }
        Some(array.map(|x| unsafe { x.assume_init() }))
    }
}

pub trait IterExt: Sized + Iterator {
    fn array_chunks<const N: usize>(self) -> ArrayChunks<Self, N>;
}

impl<I: Iterator> IterExt for I {
    fn array_chunks<const N: usize>(self) -> ArrayChunks<Self, N> {
        ArrayChunks {
            iterator: self,
            remainder: None,
        }
    }
}
