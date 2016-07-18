#ifndef STACK_H
#define STACK_H
#include <iostream>       // std::cerr
#include <stdexcept>      // std::out_of_range

//*****************************************************************************
//                           DATA TYPES DEFINITION
//*****************************************************************************
struct Convex_T
{
	int		vertex_Index;
	int		edge_ID;
};

struct stack_T
{
	int				size;
	int				nElements;
	struct Convex_T	*data;
};

//*****************************************************************************
//                           STACK PUBLIC FUNCTIONS
//*****************************************************************************
int 	initialize_Stack(struct stack_T *stack, int size);
struct Convex_T		peak(struct stack_T *stack);
struct Convex_T		peak_Prev(struct stack_T *stack);
void	update_Peak(struct stack_T *stack, int edge);
void	push(struct stack_T *stack, struct Convex_T edge);
void	pop(struct stack_T *stack);
int		get_Stack_nElements(struct stack_T *stack);
void	reset_Stack(struct stack_T *stack);
void	print_Stack(struct stack_T *stack);
void	finalize_Stack(struct stack_T *stack);


//*****************************************************************************
//                           	STACK CLASS
//*****************************************************************************
template <class A_Type>
class Stack_T
{
	//***********************************************************
	// ATTRIBUTES
	//***********************************************************
	int			size;			// Max # elements in stack.
	int			nElements;		// Current # elements in stack.
	A_Type 		*data;			// Elements array.

	//***********************************************************
	// MEMBERS
	//***********************************************************
	public:

	//***********************************************************
	// Constructor/Destructor.
	//***********************************************************
	Stack_T(int s);
	~Stack_T();

	//***********************************************************
	// Interface.
	//***********************************************************
	A_Type 	peak();
	A_Type 	elementAt( int index);
	void 	update_Peak(A_Type value);
	void 	push(A_Type value);
	A_Type 	pop();
	int 	get_nElements();
	int 	get_size();
	void 	reset_Stack();
};

// Stack class definition.
// Constructor.
template <class A_Type> Stack_T<A_Type>::Stack_T(int s)
{
	// Set size and # elements in stack.
	this->size = s;
	this->nElements = 0;

	// Allocate array.
	this->data = new A_Type[s];
}

// Destructor.
template <class A_Type> Stack_T<A_Type>::~Stack_T()
{
	// Deallocate array.
	delete this->data;
}

// Return peak value.
template <class A_Type> A_Type Stack_T<A_Type>::peak()
{
	A_Type ret;			// Return value.

	try
	{
		// Get peak element.
		ret = this->data[this->nElements-1];
	}
	catch (const std::out_of_range& ex)
	{
		std::cerr << "Out of Range error: " << ex.what() << '\n';
	}

	return(ret);
}

// Return element at position index.
template <class A_Type> A_Type Stack_T<A_Type>::elementAt(int index)
{
	A_Type ret;			// Return value.

	try
	{
		// Get peak element.
		ret = this->data[this->nElements-index];
	}
	catch (const std::out_of_range& ex)
	{
		std::cerr << "Out of Range error: " << ex.what() << '\n';
	}

	return(ret);
}

// Set peak value without inserting a new value.
template <class A_Type> void Stack_T<A_Type>::update_Peak(A_Type value)
{
	try
	{
		// Set peak value.
		this->data[this->nElements-1] = value;
	}
	catch (const std::out_of_range& ex)
	{
		std::cerr << "Out of Range error: " << ex.what() << '\n';
	}
}

// Add a new element at top of the stack.
template <class A_Type> void Stack_T<A_Type>::push(A_Type value)
{
	try
	{
		// Insert a new element in the top of the stack.
		this->data[this->nElements] = value;
		this->nElements++;
	}
	catch (const std::out_of_range& ex)
	{
		std::cerr << "Out of Range error: " << ex.what() << '\n';
	}
}

// Return peak element and remove it from the stack.
template <class A_Type> A_Type Stack_T<A_Type>::pop()
{
	A_Type ret;		// Return value.

	try
	{
		// Decrease number of elements.
		this->nElements--;

		ret = this->data[this->nElements];
	}
	catch (const std::out_of_range& ex)
	{
		std::cerr << "Out of Range error: " << ex.what() << '\n';
	}

	return(ret);
}

// Return nElements value.
template <class A_Type> int Stack_T<A_Type>::get_nElements()
{
	// Return current # elements in stack.
	return(this->nElements);
}

// Return stack size.
template <class A_Type> int Stack_T<A_Type>::get_size()
{
	// Return stack size.
	return(this->size);
}

// Set stack as empty.
template <class A_Type> void Stack_T<A_Type>::reset_Stack()
{
	// Reset data of the stack (without resizing).
	this->nElements = 0;
	memset( this->data, 0, sizeof(A_Type)*this->size);
}

#endif
