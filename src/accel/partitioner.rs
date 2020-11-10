use crate::primitive::*;

pub struct PartitionInfo {
    pub index: usize,
    pub axis: u8,
}

pub trait Partitioner {
    fn new() -> Self;
    fn partition<T: Bounded3<Float>>(&self, primitives: &mut [T]) -> PartitionInfo;
}

// partitions a set of primitives on the middle point of the largest extent
pub struct MiddlePartitioner {}

impl Partitioner for MiddlePartitioner {
    fn new() -> Self {
        Self {}
    }

    fn partition<T: Bounded3<Float>>(&self, primitives: &mut [T]) -> PartitionInfo {
        let bbox = {
            let mut bbox = Bbox3f::new();
            for p in &primitives[..] {
                bbox += p.bbox();
            }
            bbox
        };
        let ((min, max), axis) = bbox.max_extent();

        PartitionInfo {
            index: partition_in_place(primitives.into_iter(), |p| {
                p.bbox().center().get(axis) < (min + max) / 2.
            }),
            axis: axis,
        }
    }
}

// TODO
pub struct SAHPartitioner {}

// stolen from https://doc.rust-lang.org/src/core/iter/traits/iterator.rs.html#97-3286
// the feature is (was) unstable
// remove this when stabilized
fn partition_in_place<'a, T: 'a>(
    ref mut it: impl DoubleEndedIterator<Item = &'a mut T>,
    ref mut predicate: impl FnMut(&T) -> bool,
) -> usize {
    #[inline]
    fn is_false<'a, T>(
        predicate: &'a mut impl FnMut(&T) -> bool,
        true_count: &'a mut usize,
    ) -> impl FnMut(&&mut T) -> bool + 'a {
        move |x| {
            let p = predicate(&**x);
            *true_count += p as usize;
            !p
        }
    }

    #[inline]
    fn is_true<T>(predicate: &mut impl FnMut(&T) -> bool) -> impl FnMut(&&mut T) -> bool + '_ {
        move |x| predicate(&**x)
    }

    // Repeatedly find the first `false` and swap it with the last `true`.
    let mut true_count = 0;
    while let Some(head) = it.find(is_false(predicate, &mut true_count)) {
        if let Some(tail) = it.rfind(is_true(predicate)) {
            std::mem::swap(head, tail);
            true_count += 1;
        } else {
            break;
        }
    }
    true_count
}

#[test]
fn test_partition_in_place() {
    let mut arr: Vec<i32> = (0..10).collect();
    let idx = partition_in_place((&mut arr).into_iter(), |&x| x % 2 == 0);
    for (i, v) in arr.iter().enumerate() {
        if i < idx {
            assert_eq!(v % 2, 0);
        } else {
            assert_eq!(v % 2, 1);
        }
    }
}
