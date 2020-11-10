use super::partitioner::*;
use crate::bbox::*;
use crate::primitive::*;
use crate::vector::*;

static MAX_PRIMS_IN_NODE: usize = 8;

pub struct BVHAccel<P: Partitioner> {
    primitives: Vec<PrimitiveHandle>,
    partitioner: P,
}

enum BVHNode {
    Internal {
        bbox: Bbox3f,
        // near child with respect to the split axis
        left_child: Box<BVHNode>,
        // far child with respect to the split axis
        right_child: Box<BVHNode>,
        split_axis: u8,
    },
    Leaf {
        bbox: Bbox3f,
        primitives_range: (usize, usize),
    },
}

impl Bounded3<Float> for BVHNode {
    fn bbox(&self) -> Bbox3f {
        match self {
            BVHNode::Internal { bbox, .. } => *bbox,
            BVHNode::Leaf { bbox, .. } => *bbox,
        }
    }
}

impl BVHNode {
    pub fn init_internal(
        left_child: Box<BVHNode>,
        right_child: Box<BVHNode>,
        split_axis: u8,
    ) -> Self {
        let bbox = (*left_child).bbox() + (*right_child).bbox();
        BVHNode::Internal {
            left_child,
            right_child,
            split_axis,
            bbox,
        }
    }

    pub fn init_leaf(primitives_range: (usize, usize), primitives: &[PrimitiveHandle]) -> Self {
        let bbox = {
            let mut bbox = Bbox3f::new();
            for p in &primitives[..] {
                bbox += (*p).bbox();
            }
            bbox
        };
        BVHNode::Leaf {
            bbox,
            primitives_range,
        }
    }
}

impl<P: Partitioner> BVHAccel<P> {
    pub fn init(primitives: Vec<PrimitiveHandle>) -> Self {
        let mut bvh_accel = Self {
            primitives,
            partitioner: P::new(),
        };
        bvh_accel.recursive_build();
        bvh_accel
    }

    fn recursive_build(&mut self) {
        let mut ordered_prims: Vec<PrimitiveHandle> = Vec::new();
        self.recursive_build_impl(0, self.primitives.len(), &mut ordered_prims);
        self.primitives = ordered_prims;
    }

    fn recursive_build_impl(
        &mut self,
        start: usize,
        end: usize,
        ordered_prims: &mut Vec<PrimitiveHandle>,
    ) -> Box<BVHNode> {
        let n_primitives = end - start;
        if n_primitives < MAX_PRIMS_IN_NODE {
            // create leaf node
            let primitives_range = (ordered_prims.len(), ordered_prims.len() + n_primitives);
            let mut bbox = Bbox3f::new();

            for p in &self.primitives[start..end] {
                ordered_prims.push(PrimitiveHandle::clone(p));
                bbox += p.bbox();
            }

            Box::new(BVHNode::init_leaf(
                primitives_range,
                &self.primitives[start..end],
            ))
        } else {
            let PartitionInfo { index, axis } =
                self.partitioner.partition(&mut self.primitives[start..end]);
            Box::new(BVHNode::init_internal(
                self.recursive_build_impl(start, index, ordered_prims),
                self.recursive_build_impl(index, end, ordered_prims),
                axis,
            ))
        }
    }
}
