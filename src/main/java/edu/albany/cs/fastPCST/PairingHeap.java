package edu.albany.cs.fastPCST;

import java.util.ArrayList;

/**
 * Pairing heap data structure
 *
 * @author baojian bzhou6@albany.edu
 */
public class PairingHeap {

    private Node root;
    private ArrayList<Node> buffer;

    //reference parameters
    public Double get_min_firstP;
    public Integer get_min_secondP;
    public Double delete_min_firstP;
    public Integer delete_min_secondP;
    public Node link_firstP;
    public Node link_secondP;

    public PairingHeap(ArrayList<Node> shared_buffer) {
        this.root = null;
        buffer = shared_buffer;
    }

    private Node link(Node node1, Node node2) {
        if (node1 == null) {
            return node2;
        }
        if (node2 == null) {
            return node1;
        }
        Node smaller_node = node2;
        Node larger_node = node1;
        if (node1.value < node2.value) {
            smaller_node = node1;
            larger_node = node2;
        }
        larger_node.sibling = smaller_node.child;
        if (larger_node.sibling != null) {
            larger_node.sibling.left_up = larger_node;
        }
        larger_node.left_up = smaller_node;
        smaller_node.child = larger_node;
        larger_node.value -= smaller_node.child_offset;
        larger_node.child_offset -= smaller_node.child_offset;
        this.link_firstP = node1;
        this.link_secondP = node2;
        return smaller_node;
    }

    public boolean is_empty() {
        return root == null;
    }

    public boolean get_min(Double value, Integer payload) {
        //there exists some situation that root is null ??? check this further ...


        if (root != null) {
            if (!root.equals(null)) {
                value = root.value;
                payload = root.payload;
                this.get_min_firstP = value;
                this.get_min_secondP = payload;
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    public Node insert(Double value, Integer payload) {
        Node new_node = new Node();
        new_node.sibling = null;
        new_node.child = null;
        new_node.left_up = null;
        new_node.value = value;
        new_node.payload = payload;
        new_node.child_offset = 0D;
        root = link(root, new_node);
        return new_node;
    }

    public void add_to_heap(Double value) {
        if (root != null) {
            root.value += value;
            root.child_offset += value;
        }
    }

    public void decrease_key(Node node, Double from_value, Double to_value) {
        Double additional_offset = from_value - node.value;
        node.child_offset += additional_offset;
        node.value = to_value;
        if (node.left_up != null) {
            if (node.left_up.child == node) {
                node.left_up.child = node.sibling;
            } else {
                node.left_up.sibling = node.sibling;
            }
            if (node.sibling != null) {
                node.sibling.left_up = node.left_up;
            }
            node.left_up = null;
            node.sibling = null;
            root = link(root, node);
            node = this.link_secondP;
        }
    }

    public Node getRoot() {
        return root;
    }

    public void print() {
        System.out.format("buffer size : %d ; the buffer info : \n", buffer.size());
        System.out.format("\n===================================");
        for (int ii = 0; ii < buffer.size(); ++ii) {
            System.out.format("value : %6f ; child_offset : %6f ; payload : %d\n", buffer.get(ii).value, buffer.get(ii).child_offset, buffer.get(ii).payload);
        }
        System.out.format("===================================\n");
    }

    public void print(Node root) {
        if (root != null) {
            System.out.format("Node [value=%f, child_offset=%f, payload=%d]\n", root.value, root.child_offset, root.payload);
            print(root.child);
            print(root.sibling);
        } else {
            return;
        }
    }

    public boolean delete_min(Double value, Integer payload) {
        if (root == null) {
            return false;
        }
        Node result = root;
        buffer = new ArrayList<Node>();
        Node cur_child = root.child;
        Node next_child = null;
        while (cur_child != null) {
            buffer.add(cur_child);
            next_child = cur_child.sibling;
            cur_child.left_up = null;
            cur_child.sibling = null;
            cur_child.value += result.child_offset;
            cur_child.child_offset += result.child_offset;
            cur_child = next_child;
        }
        int merged_children = 0;
        while (merged_children + 2 <= buffer.size()) {
            Node n = link(buffer.get(merged_children), buffer.get(merged_children + 1));
            buffer.set(merged_children, this.link_firstP);
            buffer.set(merged_children + 1, this.link_secondP);
            buffer.set(merged_children / 2, n);
            merged_children += 2;
        }
        if (merged_children != buffer.size()) {
            buffer.set(merged_children / 2, buffer.get(merged_children));
            int newSize = ((merged_children / 2) + 1) - buffer.size();
            if (buffer.size() > ((merged_children / 2) + 1)) {
                for (int i = 0; i < Math.abs(newSize); i++) {
                    buffer.remove(buffer.size() - 1);
                }
            } else {
                for (int i = 0; i < newSize; i++) {
                    buffer.add(new Node());
                }
            }
        } else {
            int newSize = (merged_children / 2) - buffer.size();
            if (buffer.size() > (merged_children / 2)) {
                for (int i = 0; i < Math.abs(newSize); i++) {
                    buffer.remove(buffer.size() - 1);
                }
            } else {
                for (int i = 0; i < newSize; i++) {
                    buffer.add(new Node());
                }
            }
        }
        if (buffer.size() > 0) {
            root = buffer.get(buffer.size() - 1);
            for (int ii = buffer.size() - 2; ii >= 0; --ii) {
                root = link(root, buffer.get(ii));
                buffer.set(ii, this.link_secondP);
            }
        } else {
            root = null;
        }
        for (int i = 0; i < buffer.size(); i++) {

        }

        this.delete_min_firstP = result.value;
        this.delete_min_secondP = result.payload;
        return true;
    }

    public PairingHeap meld(PairingHeap heap1, PairingHeap heap2) {
        PairingHeap result = new PairingHeap(buffer);
        result.root = link(heap1.root, heap2.root);
        heap1.root = null;
        heap2.root = null;
        return result;
    }

    public class Node {
        public Node sibling;
        public Node child;
        public Node left_up;
        public Double value;
        public Double child_offset;
        public Integer payload;

        @Override
        public String toString() {
            return "Node [value=" + value + ", child_offset=" + child_offset + ", payload=" + payload + "]";
        }
    }

}
