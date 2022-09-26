package project_utils;

/**
 * Double linked list, not cyclic.
 * 
 * @author kroka
 */
public class List<T> {
    ListElem head = null;
    ListElem last = null;
    private int size = 0;
    
    public List(){}
    
    public void add(T elem){
        if(head == null){
            head = new ListElem(elem);
            last = head;
        } else {
            ListElem n = new ListElem(elem);
            n.prev = last;
            last.next = n;
            last = n;
        }
        size++;
    }
    
    protected class ListElem {
        T elem;
        ListElem prev = null;
        ListElem next = null;
        ListElem(T elem){
            ListElem.this.elem = elem;
        }
    }
}
